################################################
# Refactored Clean task
#
# v1.0: 2012.10.05, U.R.V.
#
################################################

from __future__ import absolute_import

import os
import shutil
import numpy
import copy
import filecmp
import time
# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *
if is_CASA6:
    from casatasks import casalog

    from casatasks.private.imagerhelpers.imager_base import PySynthesisImager
    from casatasks.private.imagerhelpers.input_parameters import saveparams2last
    from casatasks.private.imagerhelpers.imager_parallel_continuum import PyParallelContSynthesisImager
    from casatasks.private.imagerhelpers.imager_parallel_cube import PyParallelCubeSynthesisImager
    from casatasks.private.imagerhelpers.input_parameters import ImagerParameters
    from .cleanhelper import write_tclean_history, get_func_params
    from casatools import table
    from casatools import image
    from casatools import synthesisutils
    from casatools import synthesisimager
else:
    from taskinit import *

    from imagerhelpers.imager_base import PySynthesisImager
    from imagerhelpers.imager_parallel_continuum import PyParallelContSynthesisImager
    from imagerhelpers.imager_parallel_cube import PyParallelCubeSynthesisImager
    from imagerhelpers.input_parameters import ImagerParameters
    from imagerhelpers.input_parameters import saveparams2last
    from cleanhelper import write_tclean_history, get_func_params
    table=casac.table
    synthesisimager=casac.synthesisimager
    image=casac.image
try:
    if is_CASA6:
        from casampi.MPIEnvironment import MPIEnvironment
        from casampi import MPIInterface
    else:
        from mpi4casa.MPIEnvironment import MPIEnvironment
        from mpi4casa import MPIInterface
    mpi_available = True
except ImportError:
    mpi_available = False



def imageDimensions(residname):
    """
    Given the input image, uses the ia tool to query for the number of Stokes
    planes and frequency channels.

    Inputs:
    residname       Name of the residual image on disk

    Returns:
    nstokes         Number of Stokes planes in the image
    nfreq           Number of frequency channels in the image
    """
    ia = image()

    # At this point we've already checked residname exists
    ia.open(residname)
    csys = ia.coordsys()

    shape = ia.shape()

    #  Figure out which axis is which
    stokes_axis = csys.findaxisbyname('Stokes')
    freq_axis = csys.findaxisbyname('Frequency')

    nstokes = shape[stokes_axis]
    nfreq = shape[freq_axis]

    ia.close()

    return nstokes, nfreq, stokes_axis, freq_axis


def fillSummaryMinor(residname, channo, stokes, stokes_axis, freq_axis, fullsummary):
    """
    Given the input image name, and the corresponding field, channel number, and
    Stokes plane, extract the relevant information from the image to generate a
    summaryminor dict as defined
    https://casadocs.readthedocs.io/en/stable/notebooks/synthesis_imaging.html#Minor-Cycle-Summary-Dictionary

    Inputs:
    residname       Name of the input residual image, str
    channo          Channel number to query in the image, int
    stokes          Stokes plane to query in the image, int
    stokes_axis     The axis to index for Stokes, int
    freq_axis       The axis to index for frequency, int
    fullsummary     Construct a full summary, or only a subset, bool

    Returns:
    summaryminor     Dict containing the necessary (key:value) pairs
    """

    ia = image()

    ia.open(residname)
    shape = ia.shape()

    casalog.post("Stokes axis is %d and Freq axis is %d" % (stokes_axis, freq_axis), "INFO3", "task_tclean")

    if stokes_axis == 2 and freq_axis == 3:
        blc = [0, 0, stokes, channo]
        trc = [shape[0], shape[1], stokes, channo]
    elif stokes_axis == 3 and freq_axis == 2:
        blc = [0, 0, channo, stokes]
        trc = [shape[0], shape[1], channo, stokes]

    data = ia.getchunk(blc, trc, dropdeg=True)
    mask = ia.getchunk(blc, trc, dropdeg=True, getmask=True)
    ia.close()

    peak_resid = numpy.amax(data*mask)
    if fullsummary:
        peak_resid_NM = numpy.amax(data)
        mask_sum = numpy.sum(mask)

    summaryminor = dict()
    # This entire function is only invoked in the special case of niter=0
    summaryminor['iterDone'] = [0.0,]
    summaryminor['peakRes'] = [peak_resid,]
    # model flux has to be zero because no iterations were performed
    summaryminor['modelFlux'] = [0.0,]
    # No threshold because no deconvolution done
    summaryminor['cycleThresh'] = [0.0,]

    if fullsummary:
        summaryminor['cycleStartIters'] = [0.0,]
        summaryminor['startIterDone'] = [0.0,]
        summaryminor['startPeakRes'] = [peak_resid,]
        summaryminor['startModelFlux'] = [0.0,]
        summaryminor['startPeakResNM'] = [peak_resid_NM,]
        summaryminor['peakResNM'] = [peak_resid_NM,]
        summaryminor['masksum'] = [mask_sum,]
        summaryminor['mpiServer'] = [0.0,]
        summaryminor['stopCode'] = [3,]


    return summaryminor


def constructSummaryMinor(paramList):
    """
    Constructs and populates a nested dictionary containing the summaryMinor()
    information required by tclean.

    This is a duplicate of the existing summaryMinor() class in C++, but is
    being re-written here to avoid issues with modifying iterControl.

    Inputs:
    paramList        Object that contains all imaging parameters

    Returns:
    summaryMinor    Nested dictionary with the summaryMinor keys, dict
    """

    impars = paramList.allimpars
    decpars = paramList.alldecpars

    # Each field is stored as a different key in impars
    nfields = len(impars.keys())
    casalog.post("Number of input fields in IMPARS is %d" % nfields, "INFO3", "task_tclean")

    summaryminor = dict()
    for ff in range(nfields):
        residname=impars[str(ff)]['imagename']+'.residual.tt0' if(os.path.exists(impars[str(ff)]['imagename']+'.residual.tt0')) else impars[str(ff)]['imagename']+'.residual'
        casalog.post("Residname %s " % residname, "INFO3", "task_tclean")
        fullsummary = decpars[str(ff)]['fullsummary']
        nstokes, nfreq, stokes_axis, freq_axis = imageDimensions(residname)

        summaryminor[ff] = dict()

        for cc in range(nfreq):
            summaryminor[ff][cc] = dict()
            for ss in range(nstokes):
                summaryminor[ff][cc][ss] = fillSummaryMinor(residname, cc, ss, stokes_axis, freq_axis, fullsummary)

    return summaryminor



def constructResidualDict(paramList):
    """
    Construct the residual dictionary given the input image name and associated parameters. This residual dictionary is meant to
    duplicate that generated by imager.getSummary(fullsummary) for the special
    case of niter = 0 to avoid initializing the deconvolver.

    Inputs:
    bparm       The tclean inputs, defaults where not specified. dict

    Returns:
    retrec      The return dictionary with imaging statistics. dict
    """

    retrec = dict()

    imagename = paramList.allimpars['0']['imagename']

    residname= imagename+'.residual.tt0' if(os.path.exists(imagename +'.residual.tt0')) else imagename +'.residual'
    if not os.path.exists(residname):
        raise FileNotFoundError(f'{residname} does not exist on disk. Cannot construct tclean return dictionary.')

    # Initialize the values that don't need to inspect any images
    retrec['cleanstate'] = 'running'
    retrec['cyclefactor'] = paramList.getAllPars()['cyclefactor']
    retrec['cycleiterdone'] = 0
    retrec['cycleniter'] = 0
    retrec['cyclethreshold'] = 0
    retrec['interactiveiterdone'] = 0
    retrec['interactivemode'] = paramList.alldecpars['0']['interactive']
    retrec['interactiveniter'] = 0
    retrec['interactivethreshold'] = 0

    retrec['iterdone'] = 0
    retrec['loopgain'] = 0
    retrec['maxpsffraction'] = 0
    retrec['maxpsfsidelobe'] = 0
    retrec['minpsffraction'] = 0

    retrec['niter'] = 0
    retrec['nmajordone'] = 1
    retrec['nsigma'] = 0.0
    # stopcode 3 --> Zero iterations performed
    retrec['stopcode'] = 3

    retrec['summarymajor'] = numpy.array([0,])
    retrec['summaryminor'] = constructSummaryMinor(paramList)
    retrec['threshold'] = 0.0
    retrec['stopDescription'] = 'Zero iterations performed'


    return retrec

#if you want to save tclean.last.* from python call of tclean uncomment the decorator   
#@saveparams2last(multibackup=True) 
def tclean(
    ####### Data Selection
    vis,#='', 
    selectdata,
    field,#='', 
    spw,#='',
    timerange,#='',
    uvrange,#='',
    antenna,#='',
    scan,#='',
    observation,#='',
    intent,#='',
    datacolumn,#='corrected',

    ####### Image definition
    imagename,#='',
    imsize,#=[100,100],
    cell,#=['1.0arcsec','1.0arcsec'],
    phasecenter,#='J2000 19:59:28.500 +40.44.01.50',
    stokes,#='I',
    projection,#='SIN',
    startmodel,#='',

    ## Spectral parameters
    specmode,#='mfs',
    reffreq,#='',
    nchan,#=1,
    start,#='',
    width,#='',
    outframe,#='LSRK',
    veltype,#='',
    restfreq,#=[''],
#    sysvel,#='',
#    sysvelframe,#='',
    interpolation,#='',
    perchanweightdensity, #=''
    ## 
    ####### Gridding parameters
    gridder,#='ft',
    facets,#=1,
    psfphasecenter,#='',

    wprojplanes,#=1,

    ### PB
    vptable,
    mosweight, #=True
    aterm,#=True,
    psterm,#=True,
    wbawp ,#= True,
    conjbeams ,#= True,
    cfcache ,#= "",
    usepointing, #=false
    computepastep ,#=360.0,
    rotatepastep ,#=360.0,
    pointingoffsetsigdev ,#=[10.0],

    pblimit,#=0.01,
    normtype,#='flatnoise',

    ####### Deconvolution parameters
    deconvolver,#='hogbom',
    scales,#=[],
    nterms,#=1,
    smallscalebias,#=0.0
    fusedthreshold,#=0.0
    largestscale,#=-1

    ### restoration options
    restoration,
    restoringbeam,#=[],
    pbcor,

    ##### Outliers
    outlierfile,#='',

    ##### Weighting
    weighting,#='natural',
    robust,#=0.5,
    noise,#0.0Jy
    npixels,#=0,
#    uvtaper,#=False,
    uvtaper,#=[],


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
    fullsummary,#=False,
    nmajor,#=-1,

    ##### (new) Mask parameters
    usemask,#='user',
    mask,#='',
    pbmask,#='',
    # maskthreshold,#='',
    # maskresolution,#='',
    # nmask,#=0,

    ##### automask by multithresh
    sidelobethreshold,#=5.0,
    noisethreshold,#=3.0,
    lownoisethreshold,#=3.0,
    negativethreshold,#=0.0,
    smoothfactor,#=1.0,
    minbeamfrac,#=0.3, 
    cutthreshold,#=0.01,
    growiterations,#=100
    dogrowprune,#=True
    minpercentchange,#=0.0
    verbose, #=False
    fastnoise, #=False

    ## Misc

    restart,#=True,

    savemodel,#="none",

#    makeimages,#="auto"
    calcres,#=True,
    calcpsf,#=True,
    psfcutoff,#=0.35

    ####### State parameters
    parallel):#=False):

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

    if(deconvolver=="mtmfs" and (specmode=='cube' or specmode=='cubedata')):
        casalog.post( "The MSMFS algorithm (deconvolver='mtmfs') with specmode='cube' is not supported", "WARN", "task_tclean" )
        return


    if((specmode=='cube' or specmode=='cubedata') and (parallel==False and mpi_available and   MPIEnvironment.is_mpi_enabled) ):
        casalog.post( "Setting parameter parallel=False with specmode='cube' when launching CASA with mpi has no effect", "WARN", "task_tclean" )
        
      
    if(perchanweightdensity==False and weighting=='briggsbwtaper'):
        casalog.post( "The briggsbwtaper weighting scheme is not compatable with perchanweightdensity=False.", "WARN", "task_tclean" )
        return
        
    if((specmode=='mfs' or specmode=='cont') and weighting=='briggsbwtaper'):
        casalog.post( "The briggsbwtaper weighting scheme is not compatable with specmode='mfs' or 'cont'.", "WARN", "task_tclean" )
        return
        
    if(npixels != 0 and weighting=='briggsbwtaper'):
        casalog.post( "The briggsbwtaper weighting scheme is not compatable with npixels != 0.", "WARN", "task_tclean" )
        return


    if(facets>1 and parallel==True):
        casalog.post("Facetted imaging currently works only in serial. Please choose pure W-projection instead.","WARN","task_tclean")

    if (nmajor < -1):
        casalog.post("Negative values less than -1 for nmajor are reserved for possible future implementation", "WARN", "task_tclean")
        return

    #####################################################
    #### Construct ImagerParameters object
    #####################################################

    imager = None
    paramList = None

    # Put all parameters into dictionaries and check them.
    ##make a dictionary of parameters that ImagerParameters take

    if is_python3:
        defparm=dict(list(zip(ImagerParameters.__init__.__code__.co_varnames[1:], ImagerParameters.__init__.__defaults__)))
    else:
        defparm=dict(zip(ImagerParameters.__init__.__func__.__code__.co_varnames[1:], ImagerParameters.__init__.func_defaults))
        
    ###assign values to the ones passed to tclean and if not defined yet in tclean...
    ###assign them the default value of the constructor
    bparm={k:  inpparams[k] if k in inpparams else defparm[k]  for k in defparm.keys()}

    ###default mosweight=True is tripping other gridders as they are not
    ###expecting it to be true
    if(bparm['mosweight']==True and bparm['gridder'].find("mosaic") == -1):
        bparm['mosweight']=False

    if specmode=='mfs':
        bparm['perchanweightdensity'] = False
    
    # deprecation message
    if usemask=='auto-thresh' or usemask=='auto-thresh2':
        casalog.post(usemask+" is deprecated, will be removed in CASA 5.4.  It is recommended to use auto-multithresh instead", "WARN") 

    #paramList.printParameters()
    
    if len(pointingoffsetsigdev)>0 and pointingoffsetsigdev[0]!=0.0 and usepointing==True and gridder.count('awproj')>1:
        casalog.post("pointingoffsetsigdev will be used for pointing corrections with AWProjection", "WARN") 
#    elif usepointing==True and pointingoffsetsigdev[0] == 0:
#        casalog.post("pointingoffsetsigdev is set to zero which is an unphysical value, will proceed with the native sky pixel resolution instead". "WARN")


    ##pcube may still need to be set to True for some combination of ftmachine etc...
    #=========================================================
    concattype=''
    pcube=False
    if parallel==True and specmode!='mfs':
        if specmode!='mfs':
            pcube=False
            parallel=False
        else:
            pcube=True
    #=========================================================
    ####set the children to load c++ libraries and applicator
    ### make workers ready for c++ based mpicommands
    cppparallel=False
    if mpi_available and MPIEnvironment.is_mpi_enabled and specmode!='mfs' and not pcube:
        mint=MPIInterface.MPIInterface()
        cl=mint.getCluster()
        if(is_CASA6):
            cl._cluster.pgc("from casatools import synthesisimager", False)
            cl._cluster.pgc("si=synthesisimager()", False)
        else:
            cl._cluster.pgc("from casac import casac", False)
            cl._cluster.pgc("si=casac.synthesisimager()", False) 
        cl._cluster.pgc("si.initmpi()", False)
        cppparallel=True
        ###ignore chanchunk
        bparm['chanchunks']=1

    # catch non operational case (parallel cube tclean with interative=T)
    if pcube and interactive:
        casalog.post( "Interactive mode is not currently supported with parallel apwproject cube CLEANing, please restart by setting interactive=F", "WARN", "task_tclean" )
        return False
    #casalog.post('parameters {}'.format(bparm))    
    paramList=ImagerParameters(**bparm)

    ## Setup Imager objects, for different parallelization schemes.
    imagerInst=PySynthesisImager
    if parallel==False and pcube==False:
         imager = PySynthesisImager(params=paramList)
         imagerInst=PySynthesisImager
    elif parallel==True:
         imager = PyParallelContSynthesisImager(params=paramList)
         imagerInst=PySynthesisImager
    elif pcube==True:
         imager = PyParallelCubeSynthesisImager(params=paramList)
         imagerInst=PyParallelCubeSynthesisImager
         # virtualconcat type - changed from virtualmove to virtualcopy 2016-07-20
         #using ia.imageconcat now the name changed to copyvirtual 2019-08-12
         concattype='copyvirtual'
    else:
         casalog.post('Invalid parallel combination in doClean.', 'ERROR')
         return
    
    retrec={}

    try: 
    #if (1):
        #pdb.set_trace()
        ## Init major cycle elements
        t0=time.time();
        imager.initializeImagers()
    
        # Construct the CFCache for AWProject-class of FTMs.  For
        # other choices the following three calls become NoOps.
        # imager.dryGridding();
        # imager.fillCFCache();
        # imager.reloadCFCache();

        imager.initializeNormalizers()
        imager.setWeighting()
        t1=time.time();
        casalog.post("***Time for initializing imager and normalizers: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

        ## Init minor cycle elements
        if niter>0 or restoration==True:
            t0=time.time();
            imager.initializeDeconvolvers()
            t1=time.time();
            casalog.post("***Time for initializing deconvolver(s): "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

        ####now is the time to check estimated memory
        imager.estimatememory()
            
        if niter>0:
            t0=time.time();
            imager.initializeIterationControl()
            t1=time.time();
            casalog.post("***Time for initializing iteration controller: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
        
        ## Make PSF
        if calcpsf==True:
            t0=time.time();
             
            imager.makePSF()
            if((psfphasecenter != '') and ('mosaic' in gridder)):
                ###for some reason imager keeps the psf open delete it and recreate it afterwards
                imager.deleteTools()
                mytb=table()
                psfname=bparm['imagename']+'.psf.tt0' if(os.path.exists( bparm['imagename']+'.psf.tt0')) else bparm['imagename']+'.psf'
                mytb.open(psfname)
                miscinf=mytb.getkeyword('miscinfo')
                iminf=mytb.getkeyword('imageinfo')
                #casalog.post ('miscinfo {} {}'.format(miscinf, iminf))
                mytb.done()
                casalog.post("doing with different phasecenter psf")
                imager.unlockimages(0)
                psfParameters=paramList.getAllPars()
                psfParameters['phasecenter']=psfphasecenter
                psfParamList=ImagerParameters(**psfParameters)
                psfimager=imagerInst(params=psfParamList)
                psfimager.initializeImagers()
                psfimager.setWeighting()
                psfimager.makeImage('psf', psfParameters['imagename']+'.psf')
                psfimager.deleteTools()
                mytb.open(psfname, nomodify=False)
                mytb.putkeyword('imageinfo',iminf)
                mytb.putkeyword('miscinfo',miscinf)
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
                if niter>0 or restoration==True:
                    imager.initializeDeconvolvers() 
                if niter>0:
                    imager.initializeIterationControl()

            t1=time.time();
            if(specmode != 'mfs' and ('stand' in gridder)):
                casalog.post("***Time for making PSF and PB: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
            else:
                casalog.post("***Time for making PSF: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
            
            imager.makePB()

            t2=time.time();
            if(specmode=='mfs' and ('stand' in gridder)):
                casalog.post("***Time for making PB: "+"%.2f"%(t2-t1)+" sec", "INFO3", "task_tclean");

        if gridder in ['mosaic','awproject']:
            imager.checkPB()

        if niter >=0 : 

            ## Make dirty image
            if calcres==True:
                t0=time.time();
                imager.runMajorCycle(isCleanCycle=False)
                t1=time.time();
                casalog.post("***Time for major cycle (calcres=T): "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean"); 

            ## In case of no deconvolution iterations....
            if niter==0 and calcres==False:
                if savemodel != "none":
                    imager.predictModel()

            # Residual image needs to be computed for this to work
            if niter==0 and calcres==True:
                retrec = constructResidualDict(paramList)

            ## Do deconvolution and iterations
            if niter>0 :
                t0=time.time();
                isit = imager.hasConverged()
                imager.updateMask()
                #if((type(usemask)==str) and ('auto' in usemask)):  
                #    isit = imager.hasConverged()
                isit = imager.hasConverged()
               
                t1=time.time();
                casalog.post("***Time to update mask: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
                while ( not isit ):
                    t0=time.time();
                    ### sometimes after automasking it does not do anything
                    doneMinor=imager.runMinorCycle()
                    t1=time.time();
                    casalog.post("***Time for minor cycle: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");

                    t0=time.time();
                    if(doneMinor):
                        imager.runMajorCycle()
                    t1=time.time();
                    casalog.post("***Time for major cycle: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
                   
                    imager.updateMask()
                    t2=time.time()
                    casalog.post("***Time to update mask: "+"%.2f"%(t2-t1)+" sec", "INFO3", "task_tclean");
                    isit = imager.hasConverged() or (not doneMinor)
                    
                ## Get summary from iterbot
                #if type(interactive) != bool:
                retrec=imager.getSummary(fullsummary);
                
                if savemodel!='none' and (interactive==True or usemask=='auto-multithresh' or nsigma>0.0):
                    paramList.resetParameters()
                    if parallel and specmode=='mfs':
                        # For parallel mfs, also needs to reset the parameters for each node
                        imager.resetSaveModelParams(paramList)
                    imager.initializeImagers()
                    imager.predictModel()

            ## Restore images.
            if restoration==True:  
                t0=time.time();
                imager.restoreImages()
                t1=time.time();
                casalog.post("***Time for restoring images: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
                if pbcor==True:
                    t0=time.time();
                    imager.pbcorImages()
                    t1=time.time();
                    casalog.post("***Time for pb-correcting images: "+"%.2f"%(t1-t0)+" sec", "INFO3", "task_tclean");
        ######### niter >=0  end if 

    finally:
        ##close tools
        # needs to deletools before concat or lock waits for ever
        if imager != None:
             imager.deleteTools()
        if(cppparallel):
            ###release workers back to python mpi control
            si=synthesisimager()
            si.releasempi()
        if (pcube):
            casalog.post("running concatImages ...")
            casalog.post("Running virtualconcat (type=%s) of sub-cubes" % concattype,"INFO2", "task_tclean")
            imager.concatImages(type=concattype)
        # CAS-10721 
        #if niter>0 and savemodel != "none":
        #    casalog.post("Please check the casa log file for a message confirming that the model was saved after the last major cycle. If it doesn't exist, please re-run tclean with niter=0,calcres=False,calcpsf=False in order to trigger a 'predict model' step that obeys the savemodel parameter.","WARN","task_tclean")

    # Write history at the end, when hopefully all .workdirectory, .work.temp, etc. are gone
    # from disk, so they won't be picked up. They need time to disappear on NFS or slow hw.
    try:
        params = get_func_params(tclean, locals())
        write_tclean_history(imagename, 'tclean', params, casalog)
    except Exception as exc:
        casalog.post("Error updating history (logtable): {} ".format(exc),'WARN')

    return retrec

##################################################
