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
from casatools import synthesisdeconvolver
from casatasks import casalog, imregrid

from .imager_base import PySynthesisImager

_ia = _image()

#############################################

class PyMtmfsViaCubeSynthesisImager(PySynthesisImager):
    """ A subclass of PySynthesisImager, for specmode='mtmfs_via_cube'

    The idea is to do the major cycle with cube imaging, then convert the cube images
    to taylor term ".ttN" images, then do the minor cycle, then convert back to cubes.
    """

    def __init__(self,params):
        super().__init__(params)
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

    def runMinorCycle(self):
        # create .ttN taylor term images for the mtmfs deconvolver
        for immod in range(0,self.NF):
            self.cube2tt(immod)

        # run the mtmfs deconvolver
        super().runMinorCycle()

        # convert back to cube images for the cube major cycle
        for immod in range(0,self.NF):
            self.tt2cube(immod)

    def cube2tt(self, immod=0):
        """ Creates the necessary taylor term images.

        Outputs:
        pb.tt0
        residual.tt0..residual.tt(N-1), model.tt0..model.tt(N-1)
        psf.tt0..psf.tt(2*N-2)

        If incompatible images already exist with the same name, replace them. """
        decpars = getDecParsForImmod(immod)
        nterms = decpars['nterms']
        imagename = decpars['imagename']

        # create the .ttN images
        imgs = [('residual',nterms), ('model',nterms), ('psf',nterms*2-1)]
        for suffix, num_terms in imgs:
            basename = f"{imagename}.{suffix}"
            for N in range(num_terms):
                ttname = f"{imagename}.{suffix}.tt{N}"

                # remove the existing image, if any
                if os.path.exists(ttname):
                    # TODO possible optimization where all the data is set to 0 instead
                    shutil.rmtree(ttname)

                # create a new, blank image based off the template baseimage
                self.makeImage(template_img=basename, output_img=ttname)

        # convert them images!
        cubewt = f"{imagename}.sumwt"
        if not os.path.exists(cubewt):
            cubewt = ""
        for suffix, num_terms in imgs:
            basename = f"{imagename}.{suffix}"
            reffreq = self.allimpars[immod]['reffreq']
            dopsf = suffix == "psf"
            self.cube_to_taylor_sum(cubename=basename, cubewt=cubewt, mtname=basename, reffreq=reffreq, nterms=nterms, dopsf=dopsf):

        # special case: just copy pb
        basename = f"{imagename}.pb"
        ttname = f"{imagename}.pb.tt0"
        shutil.copytree(basename, ttname)

    def tt2cube(self, immod=0):
        """ Creates or updates the .model image with all new data obtained
        from the .model.ttN images.
        """
        decpars = getDecParsForImmod(immod)
        nterms = decpars['nterms']
        imagename = decpars['imagename']
        reffreq = self.allimpars[immod]['reffreq']
        
        # run the conversion
        self.taylor_model_to_cube(cubename=imagename, mtname=imagename, reffreq=reffreq, nterms=nterms)

    def makeImage(self, template_img='try.psf', output_img='try.zeros.psf'):
        # get the shape
        _ia.open(template_img)
        shape = ia.shape()
        csys = ia.coordsys()
        pixeltype = ia.pixeltype()
        casalog.post(f"pixeltype: {pixeltype} ({type(pixeltype)})")
        inpixels = ia.getregion()

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
        pixels = np.zeros(shape, dtype=dtype)

        # create the new outputmask
        ia.fromarray(output_img, csys=csys.torecord(), pixels=pixels, type=pixelprefix)
        ia.close()
        ia.done()

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
    def cube_to_taylor_sum(self, cubename='', cubewt='', chanwt='', mtname='',reffreq='1.5GHz',nterms=2,dopsf=False):
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
                  It's suggested that this have same suffix as cubename.
          reffreq: reference frequency, like for tclean
          nterms: number of taylor terms to fit the spectral index to
          dopsf: Signals that cubename represents a point source function, should be true if cubename ends with ".psf".
                 If true, then output 2*nterms-1 ttN images.
        """

        refnu = _qa.convert( _qa.quantity(reffreq) ,'Hz' )['value']

        # casalog.post("&&&&&&&&& REF FREQ : " + str(refnu))

        pix=[]

        num_terms=nterms

        if dopsf==True:
            num_terms=2*nterms-1
 
        for tt in range(0,num_terms):
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

        if shp[3] != len(cwt) or len(freqlist) != len(cwt):
            raise Exception("Nchan shape mismatch between cube and sumwt.")

        cwt = cwt * chanwt  ## Merge the weights and flags. 

        sumchanwt = np.sum(cwt)  ## This is a weight

        if sumchanwt==0:
            raise Exception("Weights are all zero ! ")
        else:

            for i in range(len(freqlist)):
                wt = (freqlist[i] - refnu)/refnu
                _ia.open(cubename)
                implane = _ia.getchunk(blc=[0,0,0,i],trc=[shp[0],shp[1],0,i])
                _ia.close()
                for tt in range(0,num_terms):
                    pix[tt] = pix[tt] + (wt**tt) * implane * cwt[i]

            for tt in range(0,num_terms):
                pix[tt] = pix[tt]/sumchanwt
#        ia.close()

        for tt in range(0,num_terms):
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

        freqlist = self.getFreqList(cubename+'.psf')
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

