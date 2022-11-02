from __future__ import absolute_import
import os
import math
import shutil
import string
import time
import re;
import copy
import numpy as np

from casatools import image as _image, table as _table
from casatools import synthesisdeconvolver, quanta
from casatasks import casalog, imregrid

from .imager_base import PySynthesisImager
from .input_parameters import ImagerParameters

_ia = _image()
_tb = _table()
_qa = quanta()

#############################################

class PyMtmfsViaCubeSynthesisImager(PySynthesisImager):
    """ A subclass of PySynthesisImager, for specmode='mtmfs_via_cube'

    The idea is to do the major cycle with cube imaging, then convert the cube images
    to taylor term ".ttN" images, then do the minor cycle, then convert back to cubes.
    """

    def __init__(self, params):
        super().__init__(params)
        if self.allimpars['0']['specmode'] != 'mtmfs_via_cube':
            raise RuntimeError(f"Can't use specmode {self.allimpars['0']['specmode']} with imager helper {self.__class__.__name__}!")

        # Update some settings:
        # - specmode to cube so that we run a cube major cycle
        # - deconvolver to hogbom so that major cycle doesn't get confused TODO is this necessary?
        for k in self.allimpars:
            self.allimpars[k]['specmode'] = 'cube'
            #self.allimpars[k]['deconvolver'] = 'hogbom'
        #for k in self.allgridpars:
            #self.allgridpars[k]['deconvolver'] = 'hogbom'
        #for k in self.allnormpars:
        #    self.allnormpars[k]['deconvolver'] = 'hogbom'

        self.fresh_images = []
        self.verify_dec_pars()

#############################################
    def verify_dec_pars(self):
        for immod in range(0,self.NF):
            pars = self.alldecpars[str(immod)]
            if pars['specmode'] != 'mtmfs_via_cube':
                raise RuntimeError(f"Creating instance of class {type(self).__name__} with the wrone specmode! Expected 'mtmfs_via_cube' but instead got '{pars['specmode']}'!")
            if pars['deconvolver'] != 'mtmfs':
                raise RuntimeError(f"specmode {pars['specmode']} requires 'mtmfs' deconvolver but instead got '{pars['deconvolver']}'!")

    def get_dec_pars_for_immod(self, immod):
        pars = self.alldecpars[str(immod)]
        pars['specmode'] = 'mfs'
        return pars

    def initializeDeconvolvers(self):
        for immod in range(0,self.NF):
             self.SDtools.append(synthesisdeconvolver())
             self.SDtools[immod].setupdeconvolution(decpars=self.get_dec_pars_for_immod(immod))

#############################################
    
    def check_psf(self, immod):
        self.cube2tt(immod, suffixes=['psf', 'sumwt'])
        return super().check_psf(immod)

    def get_image_name(self, immod, suffix, ttN=None):
        decpars = self.get_dec_pars_for_immod(immod)
        imagename = decpars['imagename']
        basename = lambda img: f"{img}.{suffix}"

        if suffix == 'model':
            if not os.path.exists(basename(imagename)):
                decpars = self.get_dec_pars_for_immod(immod)
                imagename = decpars['startmodel']

        bn = basename(imagename)
        if ttN != None:
            return f"{bn}.tt{ttN}"
        return bn

    def hasConverged(self):
        # create .ttN taylor term images for the mtmfs deconvolver
        for immod in range(0,self.NF):
            suffixes = ["residual", "psf", "sumwt", "model"]
            # only need to create the psf taylor term images once (shouldn't change after check_psf)
            suffixes.remove('psf')
            self.cube2tt(immod, suffixes=suffixes)

        return super().hasConverged()

    def runMinorCycle(self):
        # convert from cube to .ttN taylor term images for the mtmfs deconvolver
        for immod in range(0,self.NF):
            # Before minorcycle : Divide out the frequency-dependent PB, multiply by a common PB.
            inpcube = self.get_image_name(immod, "residual")
            pbcube = self.get_image_name(immod, "pb")
            cubewt = self.get_image_name(immod, "sumwt")
            pblimit = self.allnormpars[str(immod)]['pblimit']
            self.modify_with_pb(inpcube=inpcube, pbcube=pbcube, cubewt=cubewt, action='div', pblimit=pblimit, freqdep=True)
            self.modify_with_pb(inpcube=inpcube, pbcube=pbcube, cubewt=cubewt, action='mult', pblimit=pblimit, freqdep=False)

            suffixes = ["residual", "psf", "sumwt", "model"]
            # only need to create the psf taylor term images once (shouldn't change after check_psf)
            suffixes.remove('psf')
            self.cube2tt(immod, suffixes=suffixes)

        # run the mtmfs deconvolver
        ret = super().runMinorCycle()

        # convert back to cube images for the cube major cycle
        for immod in range(0,self.NF):
            self.tt2cube(immod)

            # After minorcycle : Divide out the common PB, Multiply by frequency-dependent PB.
            inpcube = self.get_image_name(immod, "model")
            pbcube = self.get_image_name(immod, "pb")
            cubewt = self.get_image_name(immod, "sumwt")
            pblimit = self.allnormpars[str(immod)]['pblimit']
            self.modify_with_pb(inpcube=inpcube, pbcube=pbcube, cubewt=cubewt, action='div', pblimit=pblimit, freqdep=False)
            self.modify_with_pb(inpcube=inpcube, pbcube=pbcube, cubewt=cubewt, action='mult', pblimit=pblimit, freqdep=True)

        return ret

    def cube2tt(self, immod=0, suffixes=None):
        """ Creates the necessary taylor term images.

        Args:
          immod: which image facet/outlier field to convert
          suffixes: list of images to convert, can include any of ["residual", "psf", "sumwt"]

        Outputs:
        pb.tt0
        residual.tt0..residual.tt(N-1), model.tt0..model.tt(N-1)
        psf.tt0..psf.tt(2*N-2)

        If incompatible images already exist with the same name, replace them. """
        if suffixes == None:
            suffixes = ["residual", "psf", "sumwt"]
        decpars = self.get_dec_pars_for_immod(immod)
        nterms = decpars['nterms']

        # determine which images are being converted
        imgs = [('residual',nterms), ('psf',nterms*2-1), ('sumwt',nterms*2-1)] #, ('model',nterms)]
        tmp_imgs = []
        for suffix, num_terms in imgs:
            if suffix not in suffixes:
                continue
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
                self.fresh_images.append(ttname)

        # convert them images!
        cubewt = self.get_image_name(immod, "sumwt")
        if not os.path.exists(cubewt):
            cubewt = ""
        for suffix, num_terms in imgs:
            basename = self.get_image_name(immod, suffix)
            reffreq = self.allimpars[str(immod)]['reffreq']
            dopsf = (suffix == "psf" or suffix == "sumwt")
            self.cube_to_taylor_sum(cubename=basename, cubewt=cubewt, mtname=basename, reffreq=reffreq, nterms=nterms, dopsf=dopsf)

        # special case: just copy pb
        basename = self.get_image_name(immod, "pb")
        ttname = self.get_image_name(immod, "pb", ttN=0)
        if os.path.exists(ttname):
            shutil.rmtree(ttname)
        shutil.copytree(basename, ttname)

    def tt2cube(self, immod=0):
        """ Creates or updates the .model image with all new data obtained
        from the .model.ttN images.
        """
        decpars = self.get_dec_pars_for_immod(immod)
        nterms = decpars['nterms']
        imagename = decpars['imagename']
        reffreq = self.allimpars[str(immod)]['reffreq']
        
        # run the conversion
        self.taylor_model_to_cube(cubename=imagename, mtname=imagename, reffreq=reffreq, nterms=nterms)

    def copy_image(self, template_img='try.psf', output_img='try.zeros.psf'):
        # get the shape
        _ia.open(template_img)
        shape = _ia.shape()
        csys = _ia.coordsys()
        pixeltype = _ia.pixeltype()
        _ia.close()
        _ia.done()

        # get the data type
        dtype = np.single
        pixelprefix = pixeltype[0] # for 'f'loat, 'd'ouble, or 'c'omplex
        if pixeltype == 'double':
            dtype = np.double
        elif pixeltype == 'complex':
            dtype = np.csingle
        elif pixeltype == 'dcomplex':
            dtype = np.cdouble
            pixelprefix = 'cd'

        # populate some pixels
        shape[3] = 1 # taylor term images don't use channels
        pixels = np.zeros(shape, dtype=dtype)

        # create the new outputmask
        _ia.fromarray(output_img, csys=csys.torecord(), pixels=pixels, type=pixelprefix)
        _ia.close()
        _ia.done()

################################################
    def copy_nonexistant_keywords(self, template_img='try.psf', output_img='try.zeros.psf'):
        _tb.open(template_img)
        new_kws = _tb.getkeywords()
        _tb.close()

        _tb.open(output_img, nomodify=False)
        old_kws = _tb.getkeywords()
        for kw in old_kws:
            del new_kws[kw]
        casalog.post(f"new keywords: {new_kws}\n\n\n", "SEVERE")
        _tb.putkeywords(new_kws)
        _tb.close()

################################################
    def get_freq_list(self,imname=''):
        """ Get the list of frequencies for the given image, one for each channel.

        Returns:
          list[float] The frequencies for each channel in the image, in Hz.

        From:
          sdint_helper.py
        """

        _ia.open(imname)
        csys =_ia.coordsys()
        shp = _ia.shape()
        _ia.close()

        if(csys.axiscoordinatetypes()[3] == 'Spectral'):
             restfreq = csys.referencevalue()['numeric'][3]#/1.0e+09; # convert more generally..
             freqincrement = csys.increment()['numeric'][3]# /1.0e+09;
             freqlist = [];
             for chan in range(0,shp[3]):
                   freqlist.append(restfreq + chan * freqincrement);
        elif(csys.axiscoordinatetypes()[3] == 'Tabular'):
             freqlist = (csys.torecord()['tabular2']['worldvalues']) # /1.0e+09;
        else:
             casalog.post('Unknown frequency axis. Exiting.','SEVERE');
             return False;

        csys.done()
        return freqlist

################################################
    def cube_to_taylor_sum(self, cubename='', cubewt='', chanwt=None, mtname='',reffreq='1.5GHz',nterms=2,dopsf=False):
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
        if dopsf==True:
            nterms=2*nterms-1
 
        pix=[]
        for tt in range(0,nterms):
            _ia.open(mtname+'.tt'+str(tt))
            pix.append( _ia.getchunk() )
            _ia.close()
            pix[tt].fill(0.0)

        _ia.open(cubename)
        shp = _ia.shape()
        _ia.close()

        _ia.open(cubewt)
        cwt = _ia.getchunk()[0,0,0,:]
        _ia.close()

        freqlist = self.get_freq_list(cubename)
        if reffreq == '':
            # from task_sdintimaging.py
            reffreq =str( ( freqlist[0] + freqlist[ len(freqlist)-1 ] )/2.0 ) + 'Hz'
        refnu = _qa.convert( _qa.quantity(reffreq) ,'Hz' )['value']

        if shp[3] != len(cwt) or len(freqlist) != len(cwt):
            raise Exception("Nchan shape mismatch between cube and sumwt.")

        if chanwt == None:
            chanwt = np.ones(len(freqlist), 'float')
        cwt = cwt * chanwt  ## Merge the weights and flags. 

        sumchanwt = np.sum(cwt)  ## This is a weight
        if sumchanwt==0:
            raise Exception("Weights are all zero ! ")

        for i in range(len(freqlist)):
            wt = (freqlist[i] - refnu)/refnu
            _ia.open(cubename)
            implane = _ia.getchunk(blc=[0,0,0,i],trc=[shp[0],shp[1],0,i])
            _ia.close()
            for tt in range(0,nterms):
                pix[tt] = pix[tt] + (wt**tt) * implane * cwt[i]

        for tt in range(0,nterms):
            pix[tt] = pix[tt]/sumchanwt

        for tt in range(0,nterms):
            _ia.open(mtname+'.tt'+str(tt))
            _ia.putchunk(pix[tt])
            _ia.close()

################################################
    def taylor_model_to_cube(self, cubename='', mtname='',reffreq='1.5GHz',nterms=2):
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
        if not os.path.exists(cubename+'.model'):
            shutil.copytree(cubename+'.psf', cubename+'.model')
        _ia.open(cubename+'.model')
        _ia.set(0.0)
        _ia.setrestoringbeam(remove=True)
        _ia.setbrightnessunit('Jy/pixel')
        _ia.close()

        freqlist = self.get_freq_list(cubename+'.psf')
        if reffreq == '':
            # from task_sdintimaging.py
            reffreq =str( ( freqlist[0] + freqlist[ len(freqlist)-1 ] )/2.0 ) + 'Hz'
        refnu = _qa.convert( _qa.quantity(reffreq) ,'Hz' )['value']

        pix=[]

        for tt in range(0,nterms):
            _ia.open(mtname+'.model.tt'+str(tt))
            pix.append( _ia.getchunk() )
            _ia.close()

        _ia.open(cubename+'.model')
        shp = _ia.shape()
        _ia.close()

        implane = pix[0].copy()

        for i in range(len(freqlist)):
            wt = (freqlist[i] - refnu)/refnu
            implane.fill(0.0)
            for tt in range(0,nterms):
                implane = implane + (wt**tt) * pix[tt]
            _ia.open(cubename+'.model')
            _ia.putchunk(implane, blc=[0,0,0,i])
            _ia.close()

################################################
    def modify_with_pb(self, inpcube='', pbcube='',cubewt='', chanwt=None, action='mult',pblimit=0.2, freqdep=True):
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
        casalog.post('Modify with PB : ' + action + ' with frequency dependence ' + str(freqdep), "INFO")

        freqlist = self.get_freq_list(inpcube)

        _ia.open(inpcube)
        shp=_ia.shape()
        _ia.close()

        ##############
        ### Calculate a reference Primary Beam
        ### Weighted sum of pb cube

        refchan=0
        _ia.open(pbcube)
        pbplane = _ia.getchunk(blc=[0,0,0,refchan],trc=[shp[0],shp[1],0,refchan])
        _ia.close()
        pbplane.fill(0.0)

        if freqdep==False:
            _ia.open(cubewt) # .sumwt
            cwt = _ia.getchunk()[0,0,0,:]
            _ia.close()

            if shp[3] != len(cwt) or len(freqlist) != len(cwt):
                raise Exception("Modify with PB : Nchan shape mismatch between cube and sumwt.")

            if chanwt == None:
                chanwt = np.ones(len(freqlist), 'float')
            cwt = cwt * chanwt  ## Merge the weights and flags

            sumchanwt = np.sum(cwt)

            if sumchanwt==0:
                raise Exception("Weights are all zero ! ")
                
            for i in range(len(freqlist)):
                ## Read the pb per plane
                _ia.open(pbcube)
                pbplane = pbplane + cwt[i] * _ia.getchunk(blc=[0,0,0,i],trc=[shp[0],shp[1],0,i])
                _ia.close()
                
            pbplane = pbplane / sumchanwt

        ##############


        ## Special-case for setting the PBmask to be same for all freqs
        if freqdep==False:
            shutil.copytree(pbcube, pbcube+'_tmpcopy')

        for i in range(len(freqlist)):

            ## Read the pb per plane
            if freqdep==True:
                _ia.open(pbcube)
                pbplane = _ia.getchunk(blc=[0,0,0,i],trc=[shp[0],shp[1],0,i])
                _ia.close()

            ## Make a tmp pbcube with the same pb in all planes. This is for the mask.
            if freqdep==False:
                _ia.open(pbcube+'_tmpcopy')
                _ia.putchunk(pbplane, blc=[0,0,0,i])
                _ia.close()

            _ia.open(inpcube)
            implane = _ia.getchunk(blc=[0,0,0,i],trc=[shp[0],shp[1],0,i])

            outplane = pbplane.copy()
            outplane.fill(0.0)

            if action=='mult':
                pbplane[pbplane<pblimit]=0.0
                outplane = implane * pbplane
            else:
                implane[pbplane<pblimit]=0.0
                pbplane[pbplane<pblimit]=1.0
                outplane = implane / pbplane

            _ia.putchunk(outplane, blc=[0,0,0,i])
            _ia.close()

        # if freqdep==True:
        #     ## Set a mask based on frequency-dependent PB
        #     self.add_mask(inpcube,pbcube,pblimit)
        # else:
        if freqdep==False:
            ## Set a mask based on the PB in refchan
            self.add_mask(inpcube,pbcube+'_tmpcopy',pblimit)
            shutil.rmtree(pbcube+'_tmpcopy')

################################################
    def add_mask(self, inpimage='',pbimage='',pblimit=0.2):
        """ Create a new mask called 'pbmask' and set it as a defualt mask.

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
        _ia.open(inpimage)
        defaultmaskname=_ia.maskhandler('default')[0]
        allmasknames = _ia.maskhandler('get')
        
        # casalog.post("defaultmaskname=",defaultmaskname)
        if defaultmaskname!='' and defaultmaskname!='mask0':
            _ia.calcmask(mask='"'+pbimage+'"'+'>'+str(pblimit), name=defaultmaskname);

        elif defaultmaskname=='mask0':
            if 'pbmask' in allmasknames:
                _ia.maskhandler('delete','pbmask')
            _ia.calcmask(mask='"'+pbimage+'"'+'>'+str(pblimit), name='pbmask');

        _ia.close()
        _ia.done() 

#############################################
#############################################

