import os
import shutil
import time
import copy
import numpy as np
import functools

from casatools import image as _image, table as _table
from casatools import quanta, ms
from casatasks import casalog
from casatools import synthesisutils as su
from .imager_base import PySynthesisImager
from .input_parameters import ImagerParameters
from typing import Tuple, List, Union, Optional
import pdb
_ia = _image()
_tb = _table()
_qa = quanta()
_ms = ms()
_su = su()

SW=True

#############################################
def time_func(func):
    @functools.wraps(func)
    def wrap_time(*args, **kwargs):
        t0 = time.time()
        result = func(*args, **kwargs)
        t1 = time.time()
        ####print(f'#######Function {func.__name__!r} took {(t1-t0)}s')
        return result
    return wrap_time

class PyMtmfsViaCubeSynthesisImager(PySynthesisImager):
    """A subclass of PySynthesisImager, for specmode='mvc'

    The idea is to do the major cycle with cube imaging, then convert the cube images
    to taylor term ".ttN" images, then do the minor cycle, then convert back to cubes.
    """
    @time_func
    def __init__(self, params: ImagerParameters) -> None:
        # Set up the mfs part for deconv
        mfsparams = copy.deepcopy(params)
        mfsparams.allimpars["0"]["specmode"] = "mfs"
        mfsparams.alldecpars["0"]["specmode"] = "mfs"
        # if there is a startmodel make the model
        if (
            len(mfsparams.alldecpars["0"]["startmodel"])
            == mfsparams.alldecpars["0"]["nterms"]
        ):
            self.copy_startmodel(
                mfsparams.alldecpars["0"],
                mfsparams.allimpars["0"],
                mfsparams.allnormpars["0"],
            )
            mfsparams.alldecpars["0"]["startmodel"] = ""
            mfsparams.allimpars["0"]["startmodel"] = ""
            params.alldecpars["0"]["startmodel"] = ""
            params.allimpars["0"]["startmodel"] = ""


        #################################
        super().__init__(params)
        ## print(f'self all impars {self.allimpars}, \n allvars={vars(self)}')
        if self.allimpars["0"]["specmode"] != "mvc":
            raise RuntimeError(
                f"Can't use specmode {self.allimpars['0']['specmode']} with imager helper {self.__class__.__name__}!"
            )

        ### Set up and check nchan and reffreq settings. 
        nchan = self.allimpars["0"]["nchan"]
        freqbeg, freqwidth = self.determineFreqRange()
        if nchan < 1 :
            nchan=int((freqwidth)/(0.1*freqbeg))  #gives around 10 channel for 2:1 BW
            if nchan < 5:
                nchan=mfsparams.alldecpars["0"]["nterms"] + 1
            casalog.post('Calculating nchan from the data range to be '+str(nchan),'INFO')

        ## If nchan < nterms, complain.
        in_nterms = mfsparams.alldecpars["0"]["nterms"]
        if nchan<in_nterms:
            raise RuntimeError(
                f"nchan (={nchan}) should be >= nterms ( {in_nterms} ) for valid polynomial fits to be feasible. "
            )

        if nchan>50:
            casalog.post('For mvc (mtmfs_via_cube), one usually needs only about 10 channels across the freq range, to fit Taylor polynomials of a low order','WARN')

        in_reffreq =  mfsparams.allimpars["0"]["reffreq"]     ##### NEED TO MAKE THIS WORK FOR MULTIFIELD...
        if in_reffreq == "":  ## User has not set it. Calculate default
            midfreq = freqbeg + 0.5*freqwidth  # midpoint of the freq range found by "determineFreqRange()" in Hz
            
            for k in mfsparams.allimpars:
                mfsparams.allimpars[k]["reffreq"] = str(midfreq) + "Hz"
                params.allimpars[k]["reffreq"] = str(midfreq)+"Hz"


        #print("params.impars", mfsparams.allimpars)

        rfreq = _qa.convert(mfsparams.allimpars["0"]["reffreq"] , "Hz")["value"]
        #print("REFFREQ = ",rfreq, " ----" , mfsparams.allimpars["0"]["reffreq"])
        if rfreq < freqbeg or rfreq > freqbeg+freqwidth:
            casalog.post('The reffreq of ' + mfsparams.allimpars["0"]["reffreq"]  + ' is outside the selected frequency range of ' + str(freqbeg) + ' Hz - ' + str(freqbeg+freqwidth) + ' Hz','WARN')

        self.mfsImager = PySynthesisImager(mfsparams)

            
        freqwidth = freqwidth / nchan
        #print(f"#####freqbeg={freqbeg}, freqwidth={freqwidth}, nchan={nchan} for cube")
        # Update some settings:
        # - specmode to cube so that we run a cube major cycle
        # - deconvolver to hogbom so that major cycle doesn't get confused TODO is this necessary?
        for k in self.allimpars:
            self.allimpars[k]["specmode"] = "cube"
            self.allimpars[k]["start"] = f"{freqbeg}Hz"
            self.allimpars[k]["width"] = f"{freqwidth}Hz"
            # self.alldecpars[k]['specmode']='cube'
            # this is needed for basic check...it is not used
            self.allimpars[k]["deconvolver"] = "hogbom"
        for k in self.allgridpars:
            self.allgridpars[k]["deconvolver"] = "hogbom"
            self.allgridpars[k]["interpolation"] = "nearest"
            self.allgridpars[k]["conjbeams"]=False
            self.allgridpars[k]["wbawp"]=True
        for k in self.allnormpars:
            self.allnormpars[k]["deconvolver"] = "hogbom"
        self.weightpars['usecubebriggs']=False
        self.fresh_images: List[str] = []
        self.verify_dec_pars()
        #######################################
    @time_func
    def determineFreqRange(self) -> Tuple[np.double, np.double]:
        minFreq = 1.0e20
        maxFreq = 0.0
        #pdb.set_trace()
        for msid in self.allselpars:
            msname = self.allselpars[msid]["msname"]
            spwsel = (
                self.allselpars[msid]["spw"] if (self.allselpars[msid]["spw"]) else "*"
            )
            fieldsel = (
                self.allselpars[msid]["field"]
                if (self.allselpars[msid]["field"])
                else "*"
            )
            fieldid = _ms.msseltoindex(vis=msname, spw=spwsel, field=fieldsel)["field"][0]
            _tb.open(msname)
            fieldids=_tb.getcol("FIELD_ID")
            _tb.done()
            # have to do this because advisechansel does not work for fieldids not in main
            if fieldid not in fieldids:
                fieldid=fieldids[0]
            frange = _su.advisechansel(
                msname=msname, getfreqrange=True, fieldid=fieldid, spwselection=spwsel
            )
            if minFreq > _qa.convert(frange["freqstart"], "Hz")["value"]:
                minFreq = _qa.convert(frange["freqstart"], "Hz")["value"]
            if maxFreq < _qa.convert(frange["freqend"], "Hz")["value"]:
                maxFreq = _qa.convert(frange["freqend"], "Hz")["value"]
            if(minFreq > maxFreq):
                raise Exception("Failed to determine the frequency range to build the cube from the field and spw selection")
        
        #print(f"@@@@@@@MinFreq and MaxFreq = {minFreq},      {maxFreq}")
        freqwidth = maxFreq - minFreq
        return (minFreq, freqwidth)

    #############################################
    @time_func
    def verify_dec_pars(self) -> bool:
        for immod in range(0, self.NF):
            pars = self.alldecpars[str(immod)]
            if pars["specmode"] != "mvc":
                raise RuntimeError(
                    f"Creating instance of class {type(self).__name__} with the wrong specmode! Expected 'mvc' but instead got '{pars['specmode']}'!"
                )
            if pars["deconvolver"] != "mtmfs":
                raise RuntimeError(
                    f"specmode {pars['specmode']} requires 'mtmfs' deconvolver but instead got '{pars['deconvolver']}'!"
                )
            if pars["nterms"] < 2:
                raise RuntimeError(
                    f"specmode {pars['specmode']} requires nterms >1 !"
                )
        return True
    ##############################################
    def get_dec_pars_for_immod(self, immod: int) -> dict:
        pars = self.alldecpars[str(immod)]
        # Do not do these sneaky things here ...let the user change the parameters themselves
        # pars['specmode'] = 'mfs'
        return pars

    @time_func
    def initializeNormalizers(self):
        super().initializeNormalizers()
        self.mfsImager.initializeNormalizers()
    @time_func
    def initializeDeconvolvers(self):
        # for immod in range(0,self.NF):
        #     self.SDtools.append(synthesisdeconvolver())
        #     self.SDtools[immod].setupdeconvolution(decpars=self.get_dec_pars_for_immod(immod))
        # should initialize mfs deconvolvers
        self.mfsImager.initializeDeconvolvers()

    #############################################
    @time_func
    def check_psf(self, immod):
        self.cube2tt(immod, suffixes=["psf", "sumwt"])
        return super().check_psf(immod)

    ##################################################
    @time_func
    def copy_startmodel(self, decpars: dict, impars: dict, normpars: dict) -> None:
        """
        As tclean provides capacity for startmodel to be for field 0 only we don't
        need to deal with outlier fields
        """
        # decpars = self.get_dec_pars_for_immod(0)
        # print(f"DECPARS {decpars}")
        imagename = decpars["imagename"]
        basemod = imagename + ".model"
        if len(decpars["startmodel"]) == decpars["nterms"]:
            if not os.path.exists(basemod + ".tt0"):
                refImage = imagename + ".psf.tt0"
                if not os.path.exists(refImage):
                    self.copy_image(imagename + ".psf", refImage)
                for k in range(decpars["nterms"]):
                    self.regrid_image(
                        refImage, f"{basemod}.tt{k}", decpars["startmodel"][k]
                    )
        # As this is being called before __init__ temporarily assign alldecpars
        self.alldecpars = {"0": decpars}
        self.allimpars = {"0": impars}
        self.tt2cube(0)
        inpcube = self.get_image_name(0, "model")
        pbcube = self.get_image_name(0, "pb")
        pblimit = normpars["pblimit"]
        if SW==False:
            self.modify_cubemodel_with_pb(
                modcube=inpcube, pbcube=pbcube, pbtt0=pbcube + ".tt0", pblimit=np.fabs(pblimit)
            )
        else:
            imname = self.get_dec_pars_for_immod(0)['imagename']
            t0 = time.time()
            _su.apply_freq_dep_pb(cubename=imname,mtname=imname,pblimit=np.fabs(pblimit))
            t1 = time.time()
            #print(f'#######---- Function apply_freq_dep_pb  took {(t1-t0)}s')

        del self.alldecpars
        del self.allimpars

    ##################################################
    
    def get_image_name(self, immod: int, suffix: str, ttN: Optional[int] = None):
        decpars = self.get_dec_pars_for_immod(immod)
        imagename = decpars["imagename"]
        basename = lambda img: f"{img}.{suffix}"

        bn = basename(imagename)
        if ttN is not None:
            return f"{bn}.tt{ttN}"
        return bn

    def hasConverged(self) -> bool:
        # create .ttN taylor term images for the mtmfs deconvolver
        # for immod in range(0,self.NF):
        #    suffixes = ["residual",  "sumwt"]
        # only need to create the psf taylor term images once (shouldn't change after check_psf)
        #    self.cube2tt(immod, suffixes=suffixes)
        return self.mfsImager.hasConverged()

    #############################################

    def initializeIterationControl(self) -> None:
        self.mfsImager.initializeIterationControl()

    #############################################
    @time_func
    def runMajorCycle(self, isCleanCycle: bool = True) -> None:
        # hopefully this carries all info for writing last model
        self.IBtool = self.mfsImager.IBtool
        #time0 = time.time()
        super().runMajorCycle(isCleanCycle)
        time1 = time.time()
        suffixes = ["residual"] #, "sumwt"]
        for immod in range(0, self.NF):
            inpcube = self.get_image_name(immod, "residual")
            pbcube = self.get_image_name(immod, "pb")
            cubewt = self.get_image_name(immod, "sumwt")
            pblimit = self.allnormpars[str(immod)]["pblimit"]

            ##print("USING PB2TTPB from runMajor")
            ##self.cubePB2ttPB(pbcube, pbcube + ".tt0", cubewt, np.fabs(pblimit))
            ##self.cube2tt(immod, suffixes=["pb"])

            # self.modify_with_pb(inpcube=inpcube, pbcube=pbcube, cubewt=cubewt, action='div', pblimit=pblimit, freqdep=True)
            # self.modify_with_pb(inpcube=inpcube, pbcube=pbcube, cubewt=cubewt, action='mult', pblimit=pblimit, freqdep=False)
            if  SW==False:
                self.removePBSpectralIndex(inpcube, pbcube, pbcube + ".tt0", np.fabs(pblimit))
            else:
                imname = self.get_dec_pars_for_immod(immod)['imagename']
                t0 = time.time()
                _su.remove_freq_dep_pb(cubename=imname,mtname=imname,pblimit=np.fabs(pblimit))
                t1 = time.time()
                #print(f'#######---- Function remove_freq_dep_pb  took {(t1-t0)}s')

            self.cube2tt(immod, suffixes=suffixes)
        #time2 = time.time()
        #print(f"MAKE RESidual time, core={time1-time0} s, cube2tt={time2-time1}")

    ##############################################
    @time_func
    def makePSF(self) -> None:
        # pdb.set_trace()
        time0 = time.time()
        super().makePSFCore()
        time1 = time.time()
        #####have to ensure the pb is made
        super().makePB()
        pblimit = self.allnormpars['0']["pblimit"]
        cubewt = self.get_image_name(0, "sumwt")
        pbcube = self.get_image_name(0, "pb")

        #print("USING PB2TTPB from makePSF")
        #self.cubePB2ttPB(pbcube, pbcube + ".tt0", cubewt, np.fabs(pblimit))
        #suffixes = ["psf", "sumwt", "weight"]

        suffixes = ["pb","psf", "sumwt"]
        for immod in range(0, self.NF):
            self.cube2tt(immod, suffixes=suffixes)
            #now that we donot call dividepsfbyweight which did the fitting
            #we have to do it now explicitly
            self.mfsImager.PStools[immod].makepsfbeamset()
       
#        for immod in range(0, self.NF):
#            self.mfsImager.PStools[immod].gatherpsfweight()
#            self.mfsImager.PStools[immod].dividepsfbyweight()
        time2 = time.time()
        #print(f"MAKE psf time, core={time1-time0} s, cube2tt={time2-time1}")

    ###############################################
    @time_func
    def runMinorCycle(self) -> bool:
        # convert from cube to .ttN taylor term images for the mtmfs deconvolver
        time0 = time.time()
        # for immod in range(0,self.NF):
        # Before minorcycle : Divide out the frequency-dependent PB, multiply by a common PB.
        # inpcube = self.get_image_name(immod, "residual")
        # pbcube = self.get_image_name(immod, "pb")
        # cubewt = self.get_image_name(immod, "sumwt")
        # pblimit = self.allnormpars[str(immod)]['pblimit']
        # self.modify_with_pb(inpcube=inpcube, pbcube=pbcube, cubewt=cubewt, action='div', pblimit=pblimit, freqdep=True)
        # self.modify_with_pb(inpcube=inpcube, pbcube=pbcube, cubewt=cubewt, action='mult', pblimit=pblimit, freqdep=False)

        # suffixes = ["residual", "psf", "sumwt", "model"]
        # only need to create the psf taylor term images once (shouldn't change after check_psf)
        # suffixes.remove('psf')
        # suffixes=['model']
        # self.cube2tt(immod, suffixes=suffixes)

        # run the mtmfs deconvolver
        ret = self.mfsImager.runMinorCycle()
        time1 = time.time()
        # convert back to cube images for the cube major cycle
        for immod in range(0, self.NF):
            self.tt2cube(immod)

            # After minorcycle : Divide out the common PB, Multiply by frequency-dependent PB.
            inpcube = self.get_image_name(immod, "model")
            pbcube = self.get_image_name(immod, "pb")
            # cubewt = self.get_image_name(immod, "sumwt")
            pblimit = self.allnormpars[str(immod)]["pblimit"]
            # self.modify_with_pb(inpcube=inpcube, pbcube=pbcube, cubewt=cubewt, action='div', pblimit=pblimit, freqdep=False)
            # self.modify_with_pb(inpcube=inpcube, pbcube=pbcube, cubewt=cubewt, action='mult', pblimit=pblimit, freqdep=True)
            if SW==False:
                self.modify_cubemodel_with_pb(
                    modcube=inpcube, pbcube=pbcube, pbtt0=pbcube + ".tt0", pblimit=np.fabs(pblimit)
                )
            else:
                imname = self.get_dec_pars_for_immod(immod)['imagename']
                t0 = time.time()
                _su.apply_freq_dep_pb(cubename=imname,mtname=imname,pblimit=np.fabs(pblimit))
                t1 = time.time()
                #print(f'#######---- Function apply_freq_dep_pb  took {(t1-t0)}s')

        time2 = time.time()
        #print(f"Minorcycle time, minor={time1-time0} s, tt2cube={time2-time1}")
        return ret

    #############################################################
    @time_func
    def updateMask(self) -> None:
        return self.mfsImager.updateMask()

    ###########################################################
    def getSummary(self, fignum : int =1) -> dict:
        return self.mfsImager.getSummary(fignum)

    ##########################################################
    @time_func
    def restoreImages(self) -> None:
        self.mfsImager.restoreImages()

    ####################################################
    @time_func
    def pbcorImages(self) -> None:
        self.mfsImager.pbcorImages()

    ####################################
    @time_func
    def cube2tt(self, immod: int = 0, suffixes: Optional[List[str]] = None) -> None:
        """Creates the necessary taylor term images.

        Args:
          immod: which image facet/outlier field to convert
          suffixes: list of images to convert, can include any of ["residual", "psf", "sumwt"]

        Outputs:
        pb.tt0
        residual.tt0..residual.tt(N-1), model.tt0..model.tt(N-1)
        psf.tt0..psf.tt(2*N-2)

        If incompatible images already exist with the same name, replace them."""
        time0 = time.time()
        if suffixes is None:
            suffixes = ["pb", "residual", "psf", "sumwt"]
        decpars = self.get_dec_pars_for_immod(immod)
        nterms = decpars["nterms"]
        pblimit = self.allnormpars[str(immod)]["pblimit"]
        # determine which images are being converted
        imgs = [
            ("pb",1),
            ("residual", nterms),
            ("psf", nterms * 2 - 1),
            ("sumwt", nterms * 2 - 1)
            #            ("weight", 1),
        ]  # , ('model',nterms)]
        tmp_imgs = []
        for suffix, num_terms in imgs:
            if suffix not in suffixes:
                continue
            if os.path.exists(self.get_image_name(immod, suffix)):
                tmp_imgs.append((suffix, num_terms))
        imgs = tmp_imgs

        # create the .ttN images
        for suffix, num_terms in imgs:
            basename = self.get_image_name(immod, suffix)
            for N in range(num_terms):
                ttname = self.get_image_name(immod, suffix, ttN=N)

                # remove the existing image, if any
                if os.path.exists(ttname):
                    # only create the images once per execution
                    if ttname in self.fresh_images:
                        continue
                    # TODO possible optimization where all the data is set to 0 instead
                    shutil.rmtree(ttname)

                # create a new, blank image based off the template baseimage
                self.copy_image(template_img=basename, output_img=ttname)
                self.copy_nonexistant_keywords(template_img=basename, output_img=ttname)
                dopsf = suffix == "psf" or suffix == "sumwt"
                if not dopsf and pblimit>0.0:
                    self.add_mask(
                        ttname, self.get_image_name(immod, "pb", ttN=0), np.fabs(pblimit)
                    )
                self.fresh_images.append(ttname)

        # convert them images!
        cubewt = self.get_image_name(immod, "sumwt")
        chanweight = None
        if not os.path.exists(cubewt):
            cubewt = ""
        else:
            _ia.open(cubewt)
            # nchan = _ia.shape()[3]
            chanweight = _ia.getchunk()[0, 0, 0, :]
            _ia.done()
        # print(f'IMGS={imgs}')
        for suffix, num_terms in imgs:
            basename = self.get_image_name(immod, suffix)
            reffreq = self.allimpars[str(immod)]["reffreq"]
            dopsf = suffix == "psf" or suffix == "sumwt"
            chwgt = None if suffix != "weight" else chanweight
            if SW==False:
                self.cube_to_taylor_sum(
                    cubename=basename,
                    cubewt=cubewt,
                    chanwt=chwgt,
                    mtname=basename,
                    reffreq=reffreq,
                    nterms=num_terms,
                    dopsf=dopsf,
                )
            else:
                imname = self.get_dec_pars_for_immod(immod)['imagename']
                if suffix == "psf":
                    imtype=0   ## num_terms should be 2*nterms-1
                if suffix == "residual":
                    imtype=1   ## num_terms should be nterms
                if suffix == 'pb':   ## may not be used..... 
                    imtype=2   ## num_terms is 1 (for now)
                if suffix == "sumwt":
                    imtype=3
                t0=time.time()
                _su.cube_to_taylor_sum(cubename=imname,mtname=imname,nterms=nterms,reffreq=reffreq,imtype=imtype,pblimit=np.fabs(pblimit))
                #print(f"Imname : {imname} , suffix : {suffix} and type : {imtype} with pblimit : {pblimit}")
                t1 = time.time()
                #print(f'#######---- Function cube_to_taylor_sum  took {(t1-t0)}s')

                
            if not dopsf and pblimit>0.0:
                for theTerm in range(num_terms):
                    ttname = self.get_image_name(immod, suffix, ttN=theTerm)
                    # print(f'Adding masks to {ttname}')
                    self.add_mask(
                        ttname, self.get_image_name(immod, "pb", ttN=0), np.fabs(pblimit)
                    )
        time1 = time.time()
        #print(f"Time taken in cube2tt={time1-time0}")
        # special case: just copy pb
        # basename = self.get_image_name(immod, "pb")
        # ttname = self.get_image_name(immod, "pb", ttN=0)
        # if ((not os.path.exists(ttname)) and os.path.exists(basename)):
        #    self.copy_image(template_img=basename, output_img=ttname)
        #    self.copy_nonexistant_keywords(template_img=basename, output_img=ttname)
        #    #shutil.rmtree(ttname)
        #    #shutil.copytree(basename, ttname)
        #    self.cube_to_taylor_sum(cubename=basename, cubewt=cubewt, mtname=basename, reffreq=reffreq, nterms=1, dopsf=False)
        #    _ia.open(ttname)
        #    _ia.calc(pixels="'"+ttname+"'/max('"+ttname+"')")
        #    _ia.done()

    # end
    @time_func
    def tt2cube(self, immod: int =0) -> None:
        """Creates or updates the .model image with all new data obtained
        from the .model.ttN images.
        """
        time0 = time.time()
        decpars = self.get_dec_pars_for_immod(immod)
        nterms = decpars["nterms"]
        imagename = decpars["imagename"]
        reffreq = self.allimpars[str(immod)]["reffreq"]

        # run the conversion
        if SW==False:
            self.taylor_model_to_cube(
                cubename=imagename, mtname=imagename, reffreq=reffreq, nterms=nterms
            )
        else:
            t0=time.time()
            _su.taylor_coeffs_to_cube(cubename=imagename,mtname=imagename,reffreq=reffreq,nterms=nterms)
            t1 = time.time()
            #print(f'#######---- Function taylor_coeffs_to_cube  took {(t1-t0)}s')

        time1 = time.time()
        #print(f"Time taken in tt2cube {time1-time0}")

    #########################################################################
    @time_func
    def regrid_image(self, template_image: str ="", output_image: str ="", input_image: str ="") -> None:
        """
        template_image provides the coordinatesystem and shape onto which the
        input_image is regridded to and named output_image.
        if output_image exists on disk it gets overwritten
        """
        _ia.open(template_image)
        csys = _ia.coordsys()
        shp = _ia.shape()
        _ia.done()
        _ia.open(input_image)
        _ia.regrid(
            outfile=output_image,
            shape=shp,
            csys=csys.torecord(),
            axes=[0, 1],
            overwrite=True,
            force=False,
        )
        _ia.done()

    ######################################################################
    @time_func
    def copy_image(self, template_img: str ="try.psf", output_img: str ="try.zeros.psf") -> None:
        # get the shape
        _ia.open(template_img)
        shape = _ia.shape()
        nchan = _ia.shape()[3]
        csys = _ia.coordsys()
        pixeltype = _ia.pixeltype()
        _ia.close()
        _ia.done()

        # get the data type
        dtype = np.single
        pixelprefix = pixeltype[0]  # for 'f'loat, 'd'ouble, or 'c'omplex
        if pixeltype == "double":
            dtype = np.double
        elif pixeltype == "complex":
            dtype = np.csingle
        elif pixeltype == "dcomplex":
            dtype = np.cdouble
            pixelprefix = "cd"
        # Get the frequency and BW correct
        minfreq = csys.toworld([0, 0, 0, 0])["numeric"][3]
        maxfreq = csys.toworld([0, 0, 0, nchan - 1])["numeric"][3]
        a = csys.increment(type="spectral")
        a["numeric"] *= nchan
        csys.setincrement(a, type="spectral")
        b = csys.referencevalue(type="spectral")
        b["numeric"] = (minfreq + maxfreq) / 2.0
        csys.setreferencevalue(b, type="spectral")
        csys.setreferencepixel(0.0, "spectral")
        # populate some pixels
        shape[3] = 1  # taylor term images don't use channels
        pixels = np.zeros(shape, dtype=dtype)

        # create the new outputmask
        _ia.fromarray(output_img, csys=csys.torecord(), pixels=pixels, type=pixelprefix)
        _ia.close()
        _ia.done()

    ################################################
    @time_func
    def copy_nonexistant_keywords(
        self, template_img: str ="try.psf", output_img: str ="try.zeros.psf"
    ) -> None :
        _tb.open(template_img)
        new_ii = _tb.getkeyword("imageinfo")
        _tb.close()
        if "perplanebeams" in new_ii:
            del new_ii["perplanebeams"]
        old_ii = {}
        _tb.open(output_img, nomodify=False)
        kws = _tb.keywordnames()
        if "imageinfo" in kws:
            old_ii = _tb.getkeyword("imageinfo")
        old_ii.update(new_ii)
        # for kw in old_kws:
        #    del new_kws[kw]
        # for kw in new_kws:
        #    old_kws[kw]=new_kws[kw]
        #        casalog.post(f"{template_img} to {output_img} imageinfo: {old_ii}\n\n\n", "WARN")
        _tb.putkeyword("imageinfo", old_ii)
        _tb.close()
        _ia.open(template_img)
        miscinf = _ia.miscinfo()
        _ia.done()
        _ia.open(output_img)
        newmiscinf = _ia.miscinfo()
        newmiscinf.update(miscinf)
        _ia.setmiscinfo(newmiscinf)
        _ia.done()

    ################################################
    def get_freq_list(self, imname=""):
        """Get the list of frequencies for the given image, one for each channel.

        Returns:
          list[float] The frequencies for each channel in the image, in Hz.

        From:
          sdint_helper.py
        """

        _ia.open(imname)
        csys = _ia.coordsys()
        shp = _ia.shape()
        _ia.close()

        if csys.axiscoordinatetypes()[3] == "Spectral":
            restfreq = csys.referencevalue()["numeric"][
                3
            ]  # /1.0e+09; # convert more generally..
            freqincrement = csys.increment()["numeric"][3]  # /1.0e+09;
            freqlist = []
            for chan in range(0, shp[3]):
                freqlist.append(restfreq + chan * freqincrement)
        elif csys.axiscoordinatetypes()[3] == "Tabular":
            freqlist = csys.torecord()["tabular2"]["worldvalues"]  # /1.0e+09;
        else:
            casalog.post("Unknown frequency axis. Exiting.", "SEVERE")
            return False

        csys.done()
        return freqlist

    #######################################################
    @time_func
    def cubePB2ttPB(self, cubePB="", ttPB="", sumwt="", pblimit=0.2):
        """
        convert the cube PB to an average PB
        """
        time0 = time.time()
        # special case: just copy pb
        if (not os.path.exists(ttPB)) and os.path.exists(cubePB):
            self.copy_image(template_img=cubePB, output_img=ttPB)
            self.copy_nonexistant_keywords(template_img=cubePB, output_img=ttPB)
            # shutil.rmtree(ttname)
            # shutil.copytree(basename, ttname)
            # self.cube_to_taylor_sum(cubename=basename, cubewt=cubewt, mtname=basename, reffreq=reffreq, nterms=1, dopsf=False)
            _ia.open(ttPB)
            pix = np.zeros(_ia.shape(), dtype=np.float64)
            _ia.done()
            _ia.open(cubePB)
            #print(f"STATS of cube pb {_ia.statistics()}")
            shp = _ia.shape()
            cwt = np.ones((shp[3]))
            if os.path.exists(sumwt):
                _ib = _image()
                _ib.open(sumwt)
                if _ib.shape()[2] == shp[3]:
                    cwt = _ib.getchunk(blc=[0, 0, 0, 0], trc=[0, 0, 0, shp[3] - 1])
                _ib.done()

            for k in range(shp[3]):
                # print(f'chan {k}, max {np.max(pix)}, {np.min(pix)}')
                pix += (
                    _ia.getchunk(
                        blc=[0, 0, 0, k], trc=[shp[0] - 1, shp[1] - 1, shp[2] - 1, k]
                    )
                    * cwt[k]
                )
            _ia.done()

            _ia.open(ttPB)
            _ia.putchunk(pix)
            _ia.done()
            _ia.open(ttPB)
            _ia.calc(pixels="'" + ttPB + "'/max('" + ttPB + "')")
            _ia.calcmask(mask="'" + ttPB + "' > " + str(pblimit))
            _ia.done()
            time1 = time.time()
            #print(f"Time taken in cubePB2ttPB={time1-time0}")

    ################################################
    @time_func
    def cube_to_taylor_sum(
        self,
        cubename="",
        cubewt="",
        chanwt=None,
        mtname="",
        reffreq="1.5GHz",
        nterms=2,
        dopsf=False,
    ):
        """
        Convert Cubes (output of major cycle) to Taylor weighted averages (inputs to the minor cycle)
        Input : Cube image <cubename>, with channels weighted by image <cubewt>
        Output : Set of images : <mtname>.tt0, <mtname>.tt1, etc...
        Algorithm: I_ttN = sum([   I_v * ((f-ref)/ref)**N   for f in freqs   ])

        Args:
          cubename: Name of a cube image to interpret into a set of taylor term .ttN images, eg "try.residual", "joint.cube.psf".
          cubewt: Name of a .sumwt image that contains the per-channel weighting for the interferometer image.
          chanwt: List of 0s and 1s, one per channel, to effectively disable the effect of a channel on the resulting images.
          mtname: The prefix output name, to be concatenated with ".ttN" strings, eg "try_mt.residual", "joint.multiterm.psf"
                  These images should already exist by the time this function is called.
                  It's suggested that this have the same suffix as cubename.
          reffreq: reference frequency, like for tclean
          nterms: number of taylor terms to fit the spectral index to
          dopsf: Signals that cubename represents a point source function, should be true if cubename ends with ".psf" or ".sumwt".
                 If true, then output 2*nterms-1 ttN images.

        From:
          sdint_helper.py
        """
#        if dopsf is True:
#            nterms = 2 * nterms - 1

        pix = []
        for tt in range(0, nterms):
            _ia.open(mtname + ".tt" + str(tt))
            pix.append(_ia.getchunk())
            _ia.close()
            pix[tt].fill(0.0)

        _ia.open(cubename)
        shp = _ia.shape()
        _ia.close()

        _ia.open(cubewt)
        cwt = _ia.getchunk()[0, 0, 0, :]
        _ia.close()
        # This is a problem for mosaics cwt has no meaning one should use
        # the weightimage as a sensitivity weight
        # cwt_weight = copy.deepcopy(cwt)
        cwt.fill(1.0)

        ##########

        freqlist = self.get_freq_list(cubename)
        if reffreq == "":
            # from task_sdintimaging.py
            reffreq = str((freqlist[0] + freqlist[len(freqlist) - 1]) / 2.0) + "Hz"
        refnu = _qa.convert(_qa.quantity(reffreq), "Hz")["value"]
        # print(f'REFNU={refnu}')
        if shp[3] != len(cwt) or len(freqlist) != len(cwt):
            raise Exception("Nchan shape mismatch between cube and sumwt.")

        if isinstance(chanwt, type(None)):
            chanwt = np.ones(len(freqlist), "float")
        cwt = cwt * chanwt  # Merge the weights and flags.

        sumchanwt = np.sum(cwt)  # This is a weight
        if sumchanwt == 0:
            raise Exception("Weights are all zero ! ")

        for i in range(len(freqlist)):
            wt = (freqlist[i] - refnu) / refnu
            _ia.open(cubename)
            implane = _ia.getchunk(blc=[0, 0, 0, i], trc=[shp[0], shp[1], 0, i])
            _ia.close()
            for tt in range(0, nterms):
                pix[tt] = pix[tt] + (wt**tt) * implane * cwt[i]

        for tt in range(0, nterms):
            pix[tt] = pix[tt] / sumchanwt

        for tt in range(0, nterms):
            _ia.open(mtname + ".tt" + str(tt))
            _ia.putchunk(pix[tt])
            _ia.close()

    ################################################
    @time_func
    def taylor_model_to_cube(self, cubename="", mtname="", reffreq="1.5GHz", nterms=2):
        """
        Convert Taylor coefficients (output of minor cycle) to cube (input to major cycle)
        Input : Set of images with suffix : .model.tt0, .model.tt1, etc...
        Output : Cube .model image
        Algorithm: I_v = sum([   I_ttN * ((f-ref)/ref)**N   for f in freqs   ])

        Args:
          cubename: Name of a cube image, to be conconcatenated with ".model" or ".psf", eg "try"
                  This image will be updated with the data from the set of taylor term .ttN images from mtname.
                  The "<cubename>.model" image should already exist by the time this function is called, or
                  else the "<cubename>.psf" image will be copied and used in its place.
          mtname: The prefix input name, to be concatenated with ".model.ttN" strings, eg "try"
                  These images should already exist by the time this function is called.
                  It's suggested that this have same suffix as cubename.
          reffreq: reference frequency, like for tclean
          nterms: number of taylor terms to fit the spectral index to

        From:
          sdint_helper.py
        """
        if not os.path.exists(cubename + ".model"):
            shutil.copytree(cubename + ".psf", cubename + ".model")
        _ia.open(cubename + ".model")
        _ia.set(0.0)
        _ia.setrestoringbeam(remove=True)
        _ia.setbrightnessunit("Jy/pixel")
        _ia.close()

        freqlist = self.get_freq_list(cubename + ".psf")
        if reffreq == "":
            # from task_sdintimaging.py
            reffreq = str((freqlist[0] + freqlist[len(freqlist) - 1]) / 2.0) + "Hz"
        refnu = _qa.convert(_qa.quantity(reffreq), "Hz")["value"]

        # print(f'modelREFNU= {refnu}')
        pix = []

        for tt in range(0, nterms):
            _ia.open(mtname + ".model.tt" + str(tt))
            pix.append(_ia.getchunk())
            _ia.close()

        _ia.open(cubename + ".model")
        # shp = _ia.shape()
        _ia.close()

        implane = np.zeros(pix[0].shape, dtype=type(pix[0][0, 0, 0, 0]))

        for i in range(len(freqlist)):
            wt = (freqlist[i] - refnu) / refnu
            implane.fill(0.0)
            for tt in range(0, nterms):
                implane = implane + (wt**tt) * pix[tt]
            _ia.open(cubename + ".model")
            _ia.putchunk(implane, blc=[0, 0, 0, i])
            _ia.close()

    #############################################
    @time_func
    def modify_cubemodel_with_pb(self, modcube="", pbcube="", pbtt0="", pblimit=0.2):
        """
        divide channel model by the common average beam and multiply it back by the channel beam

        """
        time0 = time.time()
        _ia.open(modcube)
        shp = _ia.shape()
        nchan = shp[3]
        _ib = _image()
        _ib.open(pbtt0)
        avPB = _ib.getchunk()
        if (avPB.shape[0] != shp[0]) or (avPB.shape[1] != shp[1]):
            _ia.done()
            _ib.done()
            raise Exception(
                f"Modify model : shape of {pbtt0} is not the same as {modcube}"
            )
        _ib.done()
        _ib.open(pbcube)
        if nchan != _ib.shape()[3]:
            _ia.done()
            _ib.done()
            raise Exception(
                f"Modify model: number of channels of {modcube} is not the same as {pbcube}"
            )
        for k in range(nchan):
            fac = _ib.getchunk(blc=[0, 0, 0, k], trc=[shp[0] - 1, shp[1] - 1, 0, k])
            fac[avPB >= pblimit] /= avPB[avPB >= pblimit]
            fac[avPB < pblimit] = 0.0
            chandat = _ia.getchunk(blc=[0, 0, 0, k], trc=[shp[0] - 1, shp[1] - 1, 0, k])
            # print(f'Shapes fac={fac.shape}, chandat={chandat.shape}')
            chandat *= fac
            _ia.putchunk(chandat, blc=[0, 0, 0, k])
        _ia.done()
        _ib.done()
        time1 = time.time()
        #print(f"Time taken in modify cube model by PB is {time1-time0}")

    ##################
    @time_func
    def removePBSpectralIndex(self, cube="", pbcube="", pbtt0="", pblimit=0.2):
        """
        divide channel image by channel beam  and multiply it back by the
        common beam

        """
        time0 = time.time()
        _ia.open(cube)
        shp = _ia.shape()
        nchan = shp[3]
        _ib = _image()
        _ib.open(pbtt0)
        avPB = _ib.getchunk()
        if (avPB.shape[0] != shp[0]) or (avPB.shape[1] != shp[1]):
            _ia.done()
            _ib.done()
            raise Exception(
                f"removePBSpectralIndex : shape of {pbtt0} is not the same as {cube}"
            )
        _ib.done()
        _ib.open(pbcube)
        if nchan != _ib.shape()[3]:
            _ib.done()
            _ia.done()
            raise Exception(
                f"removePBSpectralIndex: number of channels of {cube} is not the same as {pbcube}"
            )
        for k in range(nchan):
            fac = copy.deepcopy(avPB)
            divid = _ib.getchunk(blc=[0, 0, 0, k], trc=[shp[0] - 1, shp[1] - 1, 0, k])
            # print(f'shapes fac={fac.shape}, divid={divid.shape}')
            # print(f'chan={k}, max min divid= {np.max(divid)}, {np.min(divid)}, fac, {np.max(fac)}, {np.min(fac)}')
            fac[divid > 0.0] /= divid[divid > 0.0]
            fac[avPB < pblimit] = 0.0
            chandat = _ia.getchunk(blc=[0, 0, 0, k], trc=[shp[0] - 1, shp[1] - 1, 0, k])
            chandat *= fac
            _ia.putchunk(chandat, blc=[0, 0, 0, k])
        _ia.done()
        _ib.done()
        time1 = time.time()
        #print(f"Time taken for removeSpecIndex = {time1-time0}")

    ################################################
    @time_func
    def modify_with_pb(
        self,
        inpcube="",
        pbcube="",
        cubewt="",
        chanwt=None,
        action="mult",
        pblimit=0.2,
        freqdep=True,
    ):
        """
        Multiply or divide by the PB

        Args:
          inpcube: The cube to be modified. For example: "try.int.cube.model"
          pbcube: The primary beam to multiply/divide by. For example: "try.int.cube.pb"
          cubewt: The per-channel weight of the inpcube. For example: "try.int.cube.sumwt"
          chanwt: List of 0s and 1s, one per channel, to effectively disable the effect of a channel on the resulting images.
          action: 'mult' or 'div', to multiply by the PB or divide by it.
          pblimit: For pixels less than this value in the PB, set those same pixels in the inpcube to zero.
          freqdep: True for channel by channel, False to use a freq-independent PB from the middle of the list before/after deconvolution

        From:
          sdint_helper.py
        """
        casalog.post(
            "Modify with PB : " + action + " with frequency dependence " + str(freqdep),
            "INFO",
        )

        freqlist = self.get_freq_list(inpcube)

        _ia.open(inpcube)
        shp = _ia.shape()
        _ia.close()

        ##############
        # Calculate a reference Primary Beam
        # Weighted sum of pb cube

        refchan = 0
        _ia.open(pbcube)
        pbplane = _ia.getchunk(
            blc=[0, 0, 0, refchan], trc=[shp[0] - 1, shp[1] - 1, 0, refchan]
        )
        _ia.close()
        pbplane.fill(0.0)

        if freqdep is False:
            _ia.open(cubewt)  # .sumwt
            cwt = _ia.getchunk()[0, 0, 0, :]
            _ia.close()

            if shp[3] != len(cwt) or len(freqlist) != len(cwt):
                raise Exception(
                    "Modify with PB : Nchan shape mismatch between cube and sumwt."
                )

            if chanwt is None:
                chanwt = np.ones(len(freqlist), "float")
            cwt = cwt * chanwt  # Merge the weights and flags

            sumchanwt = np.sum(cwt)

            if sumchanwt == 0:
                raise Exception("Weights are all zero ! ")

            for i in range(len(freqlist)):
                # Read the pb per plane
                _ia.open(pbcube)
                pbplane = pbplane + cwt[i] * _ia.getchunk(
                    blc=[0, 0, 0, i], trc=[shp[0] - 1, shp[1] - 1, 0, i]
                )
                _ia.close()

            pbplane = pbplane / sumchanwt

        ##############

        # Special-case for setting the PBmask to be same for all freqs
        if freqdep is False:
            shutil.copytree(pbcube, pbcube + "_tmpcopy")

        for i in range(len(freqlist)):

            # Read the pb per plane
            if freqdep is True:
                _ia.open(pbcube)
                pbplane = _ia.getchunk(
                    blc=[0, 0, 0, i], trc=[shp[0] - 1, shp[1] - 1, 0, i]
                )
                _ia.close()

            # Make a tmp pbcube with the same pb in all planes. This is for the mask.
            if freqdep is False:
                _ia.open(pbcube + "_tmpcopy")
                _ia.putchunk(pbplane, blc=[0, 0, 0, i])
                _ia.close()

            _ia.open(inpcube)
            implane = _ia.getchunk(blc=[0, 0, 0, i], trc=[shp[0] - 1, shp[1] - 1, 0, i])

            outplane = pbplane.copy()
            outplane.fill(0.0)

            if action == "mult":
                pbplane[pbplane < pblimit] = 0.0
                outplane = implane * pbplane
            else:
                implane[pbplane < pblimit] = 0.0
                pbplane[pbplane < pblimit] = 1.0
                outplane = implane / pbplane

            _ia.putchunk(outplane, blc=[0, 0, 0, i])
            _ia.close()

        # if freqdep==True:
        #     ## Set a mask based on frequency-dependent PB
        #     self.add_mask(inpcube,pbcube,pblimit)
        # else:
        if freqdep is False:
            # Set a mask based on the PB in refchan
            self.add_mask(inpcube, pbcube + "_tmpcopy", pblimit)
            shutil.rmtree(pbcube + "_tmpcopy")

    ################################################
    @time_func
    def add_mask(self, inpimage="", pbimage="", pblimit=0.2):
        """Create a new mask called 'pbmask' and set it as a defualt mask.

        Replaces the existing mask with a new mask based on the values in the pbimage
        and pblimit. The new mask name is either 'pbmask' or the name of the existing
        default mask.

        Args:
          inpimage: image to replace the mask on
          pbimage: image used to calculate the mask values, example "try.pb"
          pblimit: values greater than this in pbimage will be included in the mask

        From:
          sdint_helper.py
        """
        if not os.path.exists(pbimage):
            return
        _ia.open(inpimage)
        defaultmaskname = _ia.maskhandler("default")[0]
        # allmasknames = _ia.maskhandler('get')

        # casalog.post("defaultmaskname=",defaultmaskname)
        # if defaultmaskname!='' and defaultmaskname!='mask0':
        #    _ia.calcmask(mask='"'+pbimage+'"'+'>'+str(pblimit), name=defaultmaskname);

        # elif defaultmaskname=='mask0':
        #    if 'pbmask' in allmasknames:
        #        _ia.maskhandler('delete','pbmask')
        #    _ia.calcmask(mask='"'+pbimage+'"'+'>'+str(pblimit), name='pbmask');
        if defaultmaskname != "":
            _ia.done()
            return
        _ia.calcmask(mask='"' + pbimage + '"' + ">" + str(pblimit))
        _ia.close()
        _ia.done()


#############################################
#############################################
