################################################
# Refactored Clean task
#
# v1.0: 2012.10.05, U.R.V.
#
################################################

import platform
import os
import shutil
import numpy
import copy
import filecmp
import time
import pdb


from casatasks import casalog

from casatasks.private.imagerhelpers.imager_base import PySynthesisImager
from casatasks.private.imagerhelpers.input_parameters import saveparams2last
from casatasks.private.imagerhelpers.imager_parallel_continuum import PyParallelContSynthesisImager
from casatasks.private.imagerhelpers.imager_parallel_cube import PyParallelCubeSynthesisImager
from casatasks.private.imagerhelpers.imager_mtmfs_via_cube import  PyMtmfsViaCubeSynthesisImager
from casatasks.private.imagerhelpers.input_parameters import ImagerParameters
from casatasks.private.imagerhelpers.imager_return_dict import ImagingDict
from .cleanhelper import write_tclean_history, get_func_params
from casatools import table
from casatools import image
from casatools import synthesisutils
from casatools import synthesisimager

try:
    from casampi.MPIEnvironment import MPIEnvironment
    from casampi import MPIInterface

    mpi_available = True
except ImportError:
    mpi_available = False


# if you want to save tclean.last.* from python call of tclean uncomment the decorator
#@saveparams2last(multibackup=True)
def tclean(
    ####### Data Selection
    vis,  # ='',
    selectdata,
    field,  # ='',
    spw,  # ='',
    timerange,  # ='',
    uvrange,  # ='',
    antenna,  # ='',
    scan,  # ='',
    observation,  # ='',
    intent,  # ='',
    datacolumn,  # ='corrected',
    ####### Image definition
    imagename,  # ='',
    imsize,  # =[100,100],
    cell,  # =['1.0arcsec','1.0arcsec'],
    phasecenter,  # ='J2000 19:59:28.500 +40.44.01.50',
    stokes,  # ='I',
    projection,  # ='SIN',
    startmodel,  # ='',
    ## Spectral parameters
    specmode,  # ='mfs',
    reffreq,  # ='',
    nchan,  # =1,
    start,  # ='',
    width,  # ='',
    outframe,  # ='LSRK',
    veltype,  # ='',
    restfreq,  # =[''],
    #    sysvel,#='',
    #    sysvelframe,#='',
    interpolation,  # ='',
    perchanweightdensity,  # =''
    ##
    ####### Gridding parameters
    gridder,  # ='ft',
    facets,  # =1,
    psfphasecenter,  # ='',
    wprojplanes,  # =1,
    ### PB
    vptable,
    mosweight,  # =True
    aterm,  # =True,
    psterm,  # =True,
    wbawp,  # = True,
    conjbeams,  # = True,
    cfcache,  # = "",
    usepointing,  # =false
    computepastep,  # =360.0,
    rotatepastep,  # =360.0,
    pointingoffsetsigdev,  # =[10.0],
    pblimit,  # =0.01,
    normtype,  # ='flatnoise',
    ####### Deconvolution parameters
    deconvolver,  # ='hogbom',
    scales,  # =[],
    nterms,  # =1,
    smallscalebias,  # =0.0
    fusedthreshold,  # =0.0
    largestscale,  # =-1
    ### restoration options
    restoration,
    restoringbeam,  # =[],
    pbcor,
    ##### Outliers
    outlierfile,  # ='',
    ##### Weighting
    weighting,  # ='natural',
    robust,  # =0.5,
    noise,  # 0.0Jy
    npixels,  # =0,
    #    uvtaper,#=False,
    uvtaper,  # =[],
    ##### Iteration control
    niter,#=0,
    gain,#=0.1,
    threshold,#=0.0,
    nsigma,#=0.0
    cycleniter,#=0,
    cyclefactor,#=1.0,
    minpsffraction,#=0.1,
    maxpsffraction,#=0.8,
    interactive,#=False,
    nmajor,#=-1,
    fullsummary,#=False,

    ##### (new) Mask parameters
    usemask,  # ='user',
    mask,  # ='',
    pbmask,  # ='',
    # maskthreshold,#='',
    # maskresolution,#='',
    # nmask,#=0,
    ##### automask by multithresh
    sidelobethreshold,  # =5.0,
    noisethreshold,  # =3.0,
    lownoisethreshold,  # =3.0,
    negativethreshold,  # =0.0,
    smoothfactor,  # =1.0,
    minbeamfrac,  # =0.3,
    cutthreshold,  # =0.01,
    growiterations,  # =100
    dogrowprune,  # =True
    minpercentchange,  # =0.0
    verbose,  # =False
    fastnoise,  # =False
    ## Misc
    restart,  # =True,
    savemodel,  # ="none",
    #    makeimages,#="auto"
    calcres,  # =True,
    calcpsf,  # =True,
    psfcutoff,  # =0.35
    ####### State parameters
    parallel,
):  # =False):

    #####################################################
    #### Sanity checks and controls
    #####################################################
    ### Move these checks elsewhere ? 
    inpparams=locals().copy()
#    saveinputs(inpparams)
    ###now deal with parameters which are not the same name
    inpparams['msname']= inpparams.pop('vis')
    inpparams['timestr']= inpparams.pop('timerange')
    inpparams['uvdist']= inpparams.pop('uvrange')
    inpparams['obs']= inpparams.pop('observation')
    inpparams['state']= inpparams.pop('intent')
    inpparams['loopgain']=inpparams.pop('gain')
    inpparams['scalebias']=inpparams.pop('smallscalebias')
    #
    # Force chanchunks=1 always now (CAS-13400)
    inpparams['chanchunks']=1

    if specmode=='cont':
        specmode='mfs'
        inpparams['specmode']='mfs'
#    if specmode=='mfs' and nterms==1 and deconvolver == "mtmfs":
#        casalog.post( "The MTMFS deconvolution algorithm (deconvolver='mtmfs') needs nterms>1.Please set nterms=2 (or more). ", "WARN", "task_tclean" )
#        return

    # Force chanchunks=1 always now (CAS-13400)
    inpparams["chanchunks"] = 1

    if specmode == "cont":
        specmode = "mfs"
        inpparams["specmode"] = "mfs"
    #    if specmode=='mfs' and nterms==1 and deconvolver == "mtmfs":
    #        casalog.post( "The MTMFS deconvolution algorithm (deconvolver='mtmfs') needs nterms>1.Please set nterms=2 (or more). ", "WARN", "task_tclean" )
    #        return

    if deconvolver == "mtmfs" and (specmode == "cube" or specmode == "cubedata"):
        casalog.post(
            "The MSMFS algorithm (deconvolver='mtmfs') with specmode='cube' is not supported",
            "WARN",
            "task_tclean",
        )
        return

    if (
        (specmode == "cube" or specmode == "cubedata" or specmode == "cubesource")
    ) and (parallel == False and mpi_available and MPIEnvironment.is_mpi_enabled):
        casalog.post(
            "When CASA is launched with mpi, the parallel=False option has no effect for 'cube' imaging for gridder='mosaic','wproject','standard' and major cycles are always executed in parallel.\n",
            "WARN",
            "task_tclean",
        )
        # casalog.post( "Setting parameter parallel=False with specmode='cube' when launching CASA with mpi has no effect except for awproject.", "WARN", "task_tclean" )

    ## Part of CAS-13814, checking for the only options compatible with mtmfs_via_cube.  
    if specmode=="mvc": 
        if deconvolver != 'mtmfs' or nterms<=1 :           
            casalog.post("The specmode='mvc' option requires the deconvolver to be 'mtmfs' and 'nterms>1.",
                         "WARN",
                         "task_tclean")
            return

    ## Part of CAS-13814, moving warnings about pbcor and widebandpbcor from the C++ code to here, for better access to user-settings.
    if specmode=='mfs' and deconvolver=='mtmfs' and gridder in ['standard','mosaic'] and pbcor==True:
        casalog.post("For specmode='mfs' and deconvolver='mtmfs', the option of pbcor=True divides each restored Taylor coefficient image by the pb.tt0 image. This correction ignores the frequency-dependence of the primary beam and does not correct for PB spectral index. It is scientifically valid only for small fractional bandwidths. For more accurate wideband primary beam correction (if needed), please use one of the following options : (1) specmode='mvc' with gridder='standard' or 'mosaic' with pbcor=True,  (2) conjbeams=True and wbawp=True with gridder='awproject' and pbcor=True.",
                     "WARN",
                     "task_tclean")
    if perchanweightdensity == False and weighting == "briggsbwtaper":
        casalog.post(
            "The briggsbwtaper weighting scheme is not compatable with perchanweightdensity=False.",
            "WARN",
            "task_tclean",
        )
        
    if (specmode == "mfs" or specmode == "cont") and weighting == "briggsbwtaper":
        casalog.post(
            "The briggsbwtaper weighting scheme is not compatable with specmode='mfs' or 'cont'.",
            "WARN",
            "task_tclean",
        )
        return

    if npixels != 0 and weighting == "briggsbwtaper":
        casalog.post(
            "The briggsbwtaper weighting scheme is not compatable with npixels != 0.",
            "WARN",
            "task_tclean",
        )
        return

    if facets > 1 and parallel == True:
        casalog.post(
            "Facetted imaging currently works only in serial. Please choose pure W-projection instead.",
            "WARN",
            "task_tclean",
        )

    if nmajor < -1:
        casalog.post(
            "Negative values less than -1 for nmajor are reserved for possible future implementation",
            "WARN",
            "task_tclean",
        )
        return

    ## CAS-13814
    if (specmode == "mvc" and gridder == 'awproject' and (conjbeams==True or wbawp==False) ):
        casalog.post(
            "specmode='mvc' requires frequency-dependent primary beams to be used during cube gridding. Please set conjbeams=False and wbawp=True for the awproject gridder.",
            "WARN",
            "task_tclean",
        )
        return

        #CAS-13814
    if (specmode == "mvc" and gridder == 'mosaic' and conjbeams==True):
        casalog.post(
            "specmode='mvc' requires frequency-dependent primary beams to be used during cube gridding. Please set conjbeams=False with the mosaic gridder.",
            "WARN",
            "task_tclean",
        )
        return

    # CAS-14146
    if (specmode == "mfs" and deconvolver == 'mtmfs' and (gridder == 'mosaic' or gridder == 'awp2')):
        casalog.post(
            "Please consider using specmode=mvc with " + gridder + " gridder"
            " as this gridder does not implement conjbeams \n thus it needs a few major cycles to converge towards the correct answer",
            "WARN",
            "task_tclean",
        )

    # CAS-13581
    # XXX : Remove this once awp-hpg is released for general use
    if gridder == 'awphpg':
        casalog.post(
            "The awphpg gridder is not available for general use in the CASA 6.7.0 release. It will be made available in a future release.",
            "WARN",
            "task_tclean",
        )
 
    #####################################################
    #### Construct ImagerParameters object
    #####################################################

    imager = None
    paramList = None

    # Put all parameters into dictionaries and check them.
    ##make a dictionary of parameters that ImagerParameters take

    defparm = dict(
        list(
            zip(
                ImagerParameters.__init__.__code__.co_varnames[1:],
                ImagerParameters.__init__.__defaults__,
            )
        )
    )
    ###assign values to the ones passed to tclean and if not defined yet in tclean...
    ###assign them the default value of the constructor
    bparm = {k: inpparams[k] if k in inpparams else defparm[k] for k in defparm.keys()}

    ###default mosweight=True is tripping other gridders as they are not
    ###expecting it to be true
    if bparm["mosweight"] == True and (not bparm["gridder"] in ['mosaic', 'awp2']):
        bparm["mosweight"] = False

    if specmode == "mfs":
        bparm["perchanweightdensity"] = False

    # deprecation message
    if usemask == "auto-thresh" or usemask == "auto-thresh2":
        casalog.post(
            usemask
            + " is deprecated, will be removed in CASA 5.4.  It is recommended to use auto-multithresh instead",
            "WARN",
        )


    #paramList.printParameters()
    

    if len(pointingoffsetsigdev)>0 and pointingoffsetsigdev[0]!=0.0 and usepointing==True and gridder.count('awproj')>1:
        casalog.post("pointingoffsetsigdev will be used for pointing corrections with AWProjection", "WARN") 
#    elif usepointing==True and pointingoffsetsigdev[0] == 0:
#        casalog.post("pointingoffsetsigdev is set to zero which is an unphysical value, will proceed with the native sky pixel resolution instead". "WARN")

    if (
        len(pointingoffsetsigdev) > 0
        and pointingoffsetsigdev[0] != 0.0
        and usepointing == True
        and gridder.count("awproj") > 1
    ):
        casalog.post(
            "pointingoffsetsigdev will be used for pointing corrections with AWProjection",
            "WARN",
        )
    #    elif usepointing==True and pointingoffsetsigdev[0] == 0:
    #        casalog.post("pointingoffsetsigdev is set to zero which is an unphysical value, will proceed with the native sky pixel resolution instead". "WARN")

    ##pcube may still need to be set to True for some combination of ftmachine etc...
    # =========================================================
    concattype = ""
    pcube = False
    if parallel == True and specmode != "mfs":
        if specmode != "mfs":
            pcube = False
            parallel = False
        else:
            pcube = True
    # =========================================================
    ####set the children to load c++ libraries and applicator
    ### make workers ready for c++ based mpicommands
    cppparallel = False
    if (
        mpi_available
        and MPIEnvironment.is_mpi_enabled
        and specmode != "mfs"
        and not pcube
    ):
        mint = MPIInterface.MPIInterface()
        cl = mint.getCluster()
        cl._cluster.pgc("from casatools import synthesisimager", False)
        cl._cluster.pgc("si=synthesisimager()", False)
        cl._cluster.pgc("si.initmpi()", False)
        cppparallel = True
        ###ignore chanchunk
        bparm['chanchunks']=1
    ######awphpg case
    if(gridder=='awphpg'):
        localsi=synthesisimager()
        localsi.inithpg()
    # catch non operational case (parallel cube tclean with interative=T)
    if pcube and interactive:
        casalog.post(
            "Interactive mode is not currently supported with parallel apwproject cube CLEANing, please restart by setting interactive=F",
            "WARN",
            "task_tclean",
        )
        return False


    if interactive:
        # catch non operational case (parallel cube tclean with interative=T)
        if pcube:
            casalog.post(
                "Interactive mode is not currently supported with parallel apwproject cube CLEANing, please restart by setting interactive=F",
                "WARN",
                "task_tclean",
            )
            return False

        # Check for casaviewer, if it does not exist flag it up front for macOS
        # since casaviewer is no longer provided by default with macOS. Returning
        # False instead of throwing an exception results in:
        #
        #    RuntimeError: No active exception to reraise
        #
        # from tclean run from casashell.
        try:
            import casaviewer as __test_casaviewer
        except:
            if platform.system( ) == "Darwin":
                casalog.post(
                    "casaviewer is no longer available for macOS, for more information see: http://go.nrao.edu/casa-viewer-eol Please restart by setting interactive=F",
                    "WARN",
                    "task_tclean",
                )
                raise RuntimeError( "casaviewer is no longer available for macOS, for more information see: http://go.nrao.edu/casa-viewer-eol" )


    #casalog.post('parameters {}'.format(bparm))    
    paramList=ImagerParameters(**bparm)

    ## Setup Imager objects, for different parallelization schemes.
    imagerInst = PySynthesisImager
    if specmode == "mvc":
        imager = PyMtmfsViaCubeSynthesisImager(params=paramList)
        imagerInst = PyMtmfsViaCubeSynthesisImager
    elif parallel == False and pcube == False:
        imager = PySynthesisImager(params=paramList)
        imagerInst = PySynthesisImager
    elif parallel == True:
        imager = PyParallelContSynthesisImager(params=paramList)
        imagerInst = PySynthesisImager
    elif pcube == True:
        imager = PyParallelCubeSynthesisImager(params=paramList)
        imagerInst = PyParallelCubeSynthesisImager
        # virtualconcat type - changed from virtualmove to virtualcopy 2016-07-20
        # using ia.imageconcat now the name changed to copyvirtual 2019-08-12
        concattype = "copyvirtual"
    else:
        casalog.post("Invalid parallel combination in doClean.", "ERROR")
        return

    retrec = {}

    try:
        # if (1):
        # pdb.set_trace()
        ## Init major cycle elements
        t0 = time.time()
        imager.initializeImagers()

        # Construct the CFCache for AWProject-class of FTMs.  For
        # other choices the following three calls become NoOps.
        # imager.dryGridding();
        # imager.fillCFCache();
        # imager.reloadCFCache();

        imager.initializeNormalizers()
        imager.setWeighting()
        t1 = time.time()
        casalog.post(
            "***Time for initializing imager and normalizers: "
            + "%.2f" % (t1 - t0)
            + " sec",
            "INFO3",
            "task_tclean",
        )

        ## Init minor cycle elements
        if niter > 0 or restoration == True:
            t0 = time.time()
            imager.initializeDeconvolvers()
            t1 = time.time()
            casalog.post(
                "***Time for initializing deconvolver(s): "
                + "%.2f" % (t1 - t0)
                + " sec",
                "INFO3",
                "task_tclean",
            )

        ####now is the time to check estimated memory
        imager.estimatememory()

        if niter > 0:
            t0 = time.time()
            imager.initializeIterationControl()
            t1 = time.time()
            casalog.post(
                "***Time for initializing iteration controller: "
                + "%.2f" % (t1 - t0)
                + " sec",
                "INFO3",
                "task_tclean",
            )

        ## Make PSF

        if calcpsf==True:
            t0=time.time();
            #####TESTOO
            if(gridder=="awphpg" and specmode=="mfs"):
                imager.makePB()
            #####TESTOO 
            imager.makePSF()
            if (psfphasecenter != "") and ("mosaic" in gridder):
                ###for some reason imager keeps the psf open delete it and recreate it afterwards
                imager.deleteTools()
                mytb = table()
                psfname = (
                    bparm["imagename"] + ".psf.tt0"
                    if (os.path.exists(bparm["imagename"] + ".psf.tt0"))
                    else bparm["imagename"] + ".psf"
                )
                mytb.open(psfname)
                miscinf = mytb.getkeyword("miscinfo")
                iminf = mytb.getkeyword("imageinfo")
                # casalog.post ('miscinfo {} {}'.format(miscinf, iminf))
                mytb.done()
                casalog.post("doing with different phasecenter psf")
                imager.unlockimages(0)
                psfParameters = paramList.getAllPars()
                psfParameters["phasecenter"] = psfphasecenter
                psfParamList = ImagerParameters(**psfParameters)
                psfimager = imagerInst(params=psfParamList)
                psfimager.initializeImagers()
                psfimager.setWeighting()
                psfimager.makeImage("psf", psfParameters["imagename"] + ".psf")
                psfimager.deleteTools()
                mytb.open(psfname, nomodify=False)
                mytb.putkeyword("imageinfo", iminf)
                mytb.putkeyword("miscinfo", miscinf)
                mytb.done()
                mysu=synthesisutils()
                mysu.fitPsfBeam(imagename=bparm['imagename'],
                                nterms=(bparm['nterms']  if deconvolver=="mtmfs" else 1),
                                psfcutoff=bparm['psfcutoff'])
                imager = PySynthesisImager(params=paramList)
                imager.initializeImagers()
                imager.initializeNormalizers()
                imager.setWeighting()
                ###redo these as we destroyed things for lock issues
                ## Init minor cycle elements
                if niter > 0 or restoration == True:
                    imager.initializeDeconvolvers()
                if niter > 0:
                    imager.initializeIterationControl()

            t1 = time.time()
            if specmode != "mfs" and ("stand" in gridder):
                casalog.post(
                    "***Time for making PSF and PB: " + "%.2f" % (t1 - t0) + " sec",
                    "INFO3",
                    "task_tclean",
                )
            else:
                casalog.post(
                    "***Time for making PSF: " + "%.2f" % (t1 - t0) + " sec",
                    "INFO3",
                    "task_tclean",
                )

            imager.makePB()

            t2 = time.time()
            if specmode == "mfs" and ("stand" in gridder):
                casalog.post(
                    "***Time for making PB: " + "%.2f" % (t2 - t1) + " sec",
                    "INFO3",
                    "task_tclean",
                )

        if gridder in ["mosaic", "awproject"]:
            imager.checkPB()

        if niter >= 0:

            ## Make dirty image
            if calcres == True:
                t0 = time.time()
                imager.runMajorCycle(isCleanCycle=False)
                t1 = time.time()
                casalog.post(
                    "***Time for major cycle (calcres=T): "
                    + "%.2f" % (t1 - t0)
                    + " sec",
                    "INFO3",
                    "task_tclean",
                )

            ## In case of no deconvolution iterations....
            if niter == 0 and calcres == False:
                if savemodel != "none":
                    imager.predictModel()

            # CAS-13960 : Construct return dict for niter=0 case
            # If residual image does not exist, summaryminor will not be
            # populated.
            if niter==0:
                id = ImagingDict()
                retrec = id.construct_residual_dict(paramList)

            ## Do deconvolution and iterations
            if niter > 0:
                t0 = time.time()


                isit = imager.hasConverged()
                imager.updateMask()
                # if((type(usemask)==str) and ('auto' in usemask)):
                #    isit = imager.hasConverged()
                isit = imager.hasConverged()
                t1 = time.time()
                casalog.post(
                    "***Time to update mask: " + "%.2f" % (t1 - t0) + " sec",
                    "INFO3",
                    "task_tclean",
                )
                while not isit:
                    t0 = time.time()
                    ### sometimes after automasking it does not do anything
                    doneMinor = imager.runMinorCycle()
                    t1 = time.time()
                    casalog.post(
                        "***Time for minor cycle: " + "%.2f" % (t1 - t0) + " sec",
                        "INFO3",
                        "task_tclean",
                    )

                    t0 = time.time()
                    if doneMinor:
                        imager.runMajorCycle()
                    t1 = time.time()
                    casalog.post(
                        "***Time for major cycle: " + "%.2f" % (t1 - t0) + " sec",
                        "INFO3",
                        "task_tclean",
                    )

                    imager.updateMask()
                    t2 = time.time()
                    casalog.post(
                        "***Time to update mask: " + "%.2f" % (t2 - t1) + " sec",
                        "INFO3",
                        "task_tclean",
                    )
                    isit = imager.hasConverged() or (not doneMinor)

                ## Get summary from iterbot
                #if type(interactive) != bool:
                retrec=imager.getSummary(fullsummary);

                if savemodel!='none' and (interactive==True or usemask=='auto-multithresh' or nsigma>0.0):
                    paramList.resetParameters()
                    if parallel and specmode == "mfs":
                        # For parallel mfs, also needs to reset the parameters for each node
                        imager.resetSaveModelParams(paramList)
                    imager.initializeImagers()
                    imager.predictModel()

            ## Restore images.
            if restoration == True:
                t0 = time.time()
                imager.restoreImages()
                t1 = time.time()
                casalog.post(
                    "***Time for restoring images: " + "%.2f" % (t1 - t0) + " sec",
                    "INFO3",
                    "task_tclean",
                )
                if pbcor == True:
                    t0 = time.time()
                    imager.pbcorImages()
                    t1 = time.time()
                    casalog.post(
                        "***Time for pb-correcting images: "
                        + "%.2f" % (t1 - t0)
                        + " sec",
                        "INFO3",
                        "task_tclean",
                    )
        ######### niter >=0  end if

    finally:
        ##close tools
        # needs to deletools before concat or lock waits for ever
        if imager != None:
            imager.deleteTools()
        if cppparallel:
            ###release workers back to python mpi control
            si = synthesisimager()
            si.releasempi()
        if pcube:
            casalog.post("running concatImages ...")
            casalog.post(
                "Running virtualconcat (type=%s) of sub-cubes" % concattype,
                "INFO2",
                "task_tclean",
            )
            imager.concatImages(type=concattype)
        # CAS-10721
        # if niter>0 and savemodel != "none":
        #    casalog.post("Please check the casa log file for a message confirming that the model was saved after the last major cycle. If it doesn't exist, please re-run tclean with niter=0,calcres=False,calcpsf=False in order to trigger a 'predict model' step that obeys the savemodel parameter.","WARN","task_tclean")

    # Write history at the end, when hopefully all .workdirectory, .work.temp, etc. are gone
    # from disk, so they won't be picked up. They need time to disappear on NFS or slow hw.
    try:
        params = get_func_params(tclean, locals())
        write_tclean_history(imagename, "tclean", params, casalog)
    except Exception as exc:
        casalog.post("Error updating history (logtable): {} ".format(exc), "WARN")

    return retrec


##################################################
