from __future__ import absolute_import
import os
import math
import shutil
import string
import time
import re;
import copy
import numpy as np

from casatools import image as _image
from casatools import synthesisdeconvolver, quanta
from casatasks import casalog, imregrid

from .imager_base import PySynthesisImager
from .input_parameters import ImagerParameters

_ia = _image()
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
            self.allimpars[k]['deconvolver'] = 'hogbom'
        for k in self.allgridpars:
            self.allgridpars[k]['deconvolver'] = 'hogbom'
        for k in self.allnormpars:
            self.allnormpars[k]['deconvolver'] = 'hogbom'

        self.fresh_images = []
        self.verifyDecPars()

#############################################
    def verifyDecPars(self):
        for immod in range(0,self.NF):
            pars = self.alldecpars[str(immod)]
            if pars['specmode'] != 'mtmfs_via_cube':
                raise RuntimeError(f"Creating instance of class {type(self).__name__} with the wrone specmode! Expected 'mtmfs_via_cube' but instead got '{pars['specmode']}'!")
            if pars['deconvolver'] != 'mtmfs':
                raise RuntimeError(f"specmode {pars['specmode']} requires 'mtmfs' deconvolver but instead got '{pars['deconvolver']}'!")

    def getDecParsForImmod(self, immod):
        pars = self.alldecpars[str(immod)]
        pars['specmode'] = 'mfs'
        return pars

    def initializeDeconvolvers(self):
        for immod in range(0,self.NF):
             self.SDtools.append(synthesisdeconvolver())
             self.SDtools[immod].setupdeconvolution(decpars=self.getDecParsForImmod(immod))

#############################################
    
    def checkPSF(self, immod):
        self.cube2tt(immod, suffixes=['psf', 'sumwt'])
        return super().checkPSF(immod)

    def get_image_name(self, immod, suffix, ttN=None):
        decpars = self.getDecParsForImmod(immod)
        imagename = decpars['imagename']
        basename = f"{imagename}.{suffix}"

        if ttN != None:
            return f"{basename}.tt{ttN}"
        return basename

    def hasConverged(self):
        # create .ttN taylor term images for the mtmfs deconvolver
        for immod in range(0,self.NF):
            suffixes = ["residual", "psf", "sumwt"]
            # TODO is the model image ever used as part of hasConverged?
            # # the model doesn't exist for the first iteration
            # if not os.path.exists(self.get_image_name(immod, 'model')):
            #     suffixes.remove('model')

            # only need to create the psf taylor term images once (shouldn't change after checkPSF)
            suffixes.remove('psf')

            self.cube2tt(immod, suffixes=suffixes)

        return super().hasConverged()

    def runMinorCycle(self):
        # convert from cube to .ttN taylor term images for the mtmfs deconvolver
        for immod in range(0,self.NF):
            # Before minorcycle : Divide out the frequency-dependent PB, multiply by a common PB.
            suffixes = ["residual", "psf", "sumwt"]
            # TODO is the model image ever used as part of the minor cycle?
            # # the model doesn't exist for the first iteration
            # if not os.path.exists(self.get_image_name(immod, 'model')):
            #     suffixes.remove('model')

            # only need to create the psf taylor term images once (shouldn't change after checkPSF)
            suffixes.remove('psf')

            self.cube2tt(immod, suffixes=suffixes)

        # run the mtmfs deconvolver
        ret = super().runMinorCycle()

        # convert back to cube images for the cube major cycle
        for immod in range(0,self.NF):
            self.tt2cube(immod)

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
        decpars = self.getDecParsForImmod(immod)
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
                self.makeImage(template_img=basename, output_img=ttname)
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
        decpars = self.getDecParsForImmod(immod)
        nterms = decpars['nterms']
        imagename = decpars['imagename']
        reffreq = self.allimpars[str(immod)]['reffreq']
        
        # run the conversion
        self.taylor_model_to_cube(cubename=imagename, mtname=imagename, reffreq=reffreq, nterms=nterms)

    def makeImage(self, template_img='try.psf', output_img='try.zeros.psf'):
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
    def getFreqList(self,imname=''):
        """ Get the list of frequencies for the given image, one for each channel.

        Returns:
          list[float] The frequencies for each channel in the image, in Hz.
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

        freqlist = self.getFreqList(cubename)
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
        """
        if not os.path.exists(cubename+'.model'):
            shutil.copytree(cubename+'.psf', cubename+'.model')
        _ia.open(cubename+'.model')
        _ia.set(0.0)
        _ia.setrestoringbeam(remove=True)
        _ia.setbrightnessunit('Jy/pixel')
        _ia.close()

        freqlist = self.getFreqList(cubename+'.psf')
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

#############################################
#############################################

