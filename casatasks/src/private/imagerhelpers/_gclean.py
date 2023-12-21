########################################################################3
#  _gclean.py
#
# Copyright (C) 2021,2022,2023
# Associated Universities, Inc. Washington DC, USA.
#
# This script is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with this library; if not, write to the Free Software Foundation,
# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
#
# Correspondence concerning AIPS++ should be adressed as follows:
#        Internet email: aips2-request@nrao.edu.
#        Postal address: AIPS++ Project Office
#                        National Radio Astronomy Observatory
#                        520 Edgemont Road
#                        Charlottesville, VA 22903-2475 USA
#
import os
import asyncio
from functools import reduce
import copy
import numpy as np
import shutil
import time
import subprocess

from casatasks.private.imagerhelpers.imager_return_dict import ImagingDict
from casatasks import deconvolve, tclean

###
### import check versions
###
_GCV001 = True
_GCV002 = True
_GCV003 = True
_GCV004 = True


print("USING THIS GCLEAN")

# from casatasks.private.imagerhelpers._gclean import gclean
class gclean:
    '''gclean(...) creates a stream of convergence records which indicate
    the convergence quaility of the tclean process. The initial record
    describes the initial dirty image.
    It is designed for use with the interactive clean GUI, but it could
    be used independently. It can be used as a regular generator:
          for rec in gclean( vis='refim_point_withline.ms', imagename='test', imsize=512, cell='12.0arcsec',
                             specmode='cube', interpolation='nearest', nchan=5, start='1.0GHz', width='0.2GHz',
                             pblimit=-1e-05, deconvolver='hogbom', niter=500, cyclefactor=3, scales=[0, 3, 10] ):
              # use rec to decide when to stop, for example to check stopcode or peak residual:
              # if (rec[0] > 1) or (min(rec[1][0][0]['peakRes']) < 0.001):
              #     break
              print(rec)
    or as an async generator:
          async for rec in gclean( vis='refim_point_withline.ms', imagename='test', imsize=512, cell='12.0arcsec',
                                   specmode='cube', interpolation='nearest', nchan=5, start='1.0GHz', width='0.2GHz',
                                   pblimit=-1e-05, deconvolver='hogbom', niter=500, cyclefactor=3, scales=[0, 3, 10] ):
              # use rec to decide when to stop
              print(rec)


    See also: __next__(...) for a description of the returned rec

    TODO: do we need to preserve any hidden state between tclean calls for the iterbotsink and/or synthesisimager tools?
    '''

    def _tclean( self, *args, **kwargs ):
        """ Calls tclean records the arguments in the local history of tclean calls.

        The full tclean history for this instance can be retrieved via the cmds() method."""
        arg_s = ', '.join( map( lambda a: self._history_filter(len(self._exe_cmds), None, repr(a)), args ) )
        kw_s = ', '.join( map( lambda kv: self._history_filter(len(self._exe_cmds), kv[0], "%s=%s" % (kv[0],repr(kv[1]))), kwargs.items()) )
        if len(arg_s) > 0 and len(kw_s) > 0:
            parameters = arg_s + ", " + kw_s
        else:
            parameters = arg_s + kw_s
        self._exe_cmds.append( "tclean( %s )" % parameters )
        self._exe_cmds_per_iter[-1] += 1
        return tclean( *args, **kwargs )

    def _deconvolve( self, *args, **kwargs ):
        """ Calls deconvolve records the arguments in the local history of deconvolve calls.

        The full deconvolve history for this instance can be retrieved via the cmds() method."""
        arg_s = ', '.join( map( lambda a: self._history_filter(len(self._exe_cmds), None, repr(a)), args ) )
        kw_s = ', '.join( map( lambda kv: self._history_filter(len(self._exe_cmds), kv[0], "%s=%s" % (kv[0],repr(kv[1]))), kwargs.items()) )
        if len(arg_s) > 0 and len(kw_s) > 0:
            parameters = arg_s + ", " + kw_s
        else:
            parameters = arg_s + kw_s
        self._exe_cmds.append( "deconvolve( %s )" % parameters )
        self._exe_cmds_per_iter[-1] += 1
        return deconvolve( *args, **kwargs )

    def _remove_tree( self, directory ):
        if os.path.isdir(directory):
            shutil.rmtree(directory)
            self._exe_cmds.append( f'''shutil.rmtree( {repr(directory)} )''' )
            self._exe_cmds_per_iter[-1] += 1

    def cmds( self, history=False ):
        """ Returns the history of all tclean calls for this instance. If ``history``
        is set to True then the full history will be returned, otherwise the commands
        executed for generating the latest result are returned.
        """
        return self._exe_cmds if history else self._exe_cmds[-self._exe_cmds_per_iter[-1]:]

    def update( self, msg ):
        """ Interactive clean parameters update.

        Args:
            msg: dict with possible keys 'niter', 'cycleniter', 'nmajor', 'threshold', 'cyclefactor' and 'mask', 'niterleft','nmajorleft'
        """
        if 'niter' in msg:
            try:
                self._niter = int(msg['niter'])
            except ValueError:
                pass
        if 'cycleniter' in msg:
            try:
                self._cycleniter = int(msg['cycleniter'])
            except ValueError:
                pass
        if 'nmajor' in msg:
            try:
                self._nmajor = int(msg['nmajor'])
            except ValueError:
                pass
        #if 'niterleft' in msg:
        #    try:
        #        self._niter = int(msg['niterleft'])
        #    except ValueError:
        #        pass
        #if 'nmajorleft' in msg:
        #    try:
        #        self._nmajor = int(msg['nmajorleft'])
        #    except ValueError:
        #        pass

        if 'threshold' in msg:
            self._threshold = msg['threshold']
            self._threshold_to_float() # Convert str to float
        if 'cyclefactor' in msg:
            try:
                self._cyclefactor = int(msg['cyclefactor'])
            except ValueError:
                pass

    def _threshold_to_float(self):
        # Convert threshold from string to float if necessary
        if isinstance(self._threshold, str):
            if "mJy" in self._threshold:
                self._threshold = float(self._threshold.replace("mJy", "")) / 1e3
            elif "uJy" in self._threshold:
                self._threshold = float(self._threshold.replace("uJy", "")) / 1e6
            elif "Jy" in self._threshold:
                self._threshold = float(self._threshold.replace("Jy", ""))


    def __init__( self, vis, imagename, field='', spw='', timerange='', uvrange='', antenna='', scan='', observation='', intent='', datacolumn='corrected',
                  imsize=[100], cell=[ ], phasecenter='', stokes='I', startmodel='', specmode='cube', reffreq='', nchan=-1, start='', width='',
                  outframe='LSRK', veltype='radio', restfreq='', interpolation='linear', perchanweightdensity=True, gridder='standard', wprojplanes=int(1),
                  mosweight=True, psterm=False, wbawp=True, conjbeams=False, usepointing=False, pointingoffsetsigdev=[  ], pblimit=0.2, deconvolver='hogbom',
                  smallscalebias=0.0, niter=0, threshold='0.1Jy', nsigma=0.0, cycleniter=-1, nmajor=1, cyclefactor=1.0, minpsffraction=0.05, maxpsffraction=0.8,
                  scales=[], restoringbeam='', pbcor=False, nterms=int(2), weighting='natural', robust=float(0.5), npixels=0, gain=float(0.1),
                  sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, minbeamfrac=0.3, growiterations=75, dogrowprune=True,
                  minpercentchange=-1.0, fastnoise=True, savemodel='none', usemask='user', mask='', parallel=False,
                  history_filter=lambda index, arg, history_value: history_value ):
        self._vis = vis
        self._imagename = imagename
        self._imsize = imsize
        self._cell = cell
        self._phasecenter = phasecenter
        self._stokes = stokes
        self._startmodel = startmodel
        self._specmode = specmode
        self._reffreq = reffreq
        self._nchan = nchan
        self._start = start
        self._width = width
        self._outframe = outframe
        self._veltype = veltype
        self._restfreq = restfreq
        self._interpolation = interpolation
        self._perchanweightdensity = perchanweightdensity
        self._gridder = gridder
        self._wprojplanes = wprojplanes
        self._mosweight = mosweight
        self._psterm = psterm
        self._wbawp = wbawp
        self._conjbeams = conjbeams
        self._usepointing = usepointing
        self._pointingoffsetsigdev = pointingoffsetsigdev
        self._pblimit = pblimit
        self._deconvolver = deconvolver
        self._smallscalebias = smallscalebias
        self._niter = niter
        self._threshold = threshold
        self._cycleniter = cycleniter
        self._minpsffraction = minpsffraction
        self._maxpsffraction = maxpsffraction
        self._nsigma = nsigma
        self._nmajor = nmajor
        self._cyclefactor = cyclefactor
        self._scales = scales
        self._restoringbeam = restoringbeam
        self._pbcor = pbcor
        self._nterms = nterms
        self._exe_cmds = [ ]
        self._exe_cmds_per_iter = [ ]
        self._history_filter = history_filter
        self._finalized = False
        self._field = field
        self._spw = spw
        self._timerange = timerange
        self._uvrange = uvrange
        self._antenna = antenna
        self._scan = scan
        self._observation = observation
        self._intent = intent
        self._datacolumn = datacolumn
        self._weighting = weighting
        self._robust = robust
        self._npixels = npixels
        self._gain = gain
        self._sidelobethreshold = sidelobethreshold
        self._noisethreshold = noisethreshold
        self._lownoisethreshold = lownoisethreshold
        self._negativethreshold = negativethreshold
        self._minbeamfrac = minbeamfrac
        self._growiterations = growiterations
        self._dogrowprune = dogrowprune
        self._minpercentchange = minpercentchange
        self._fastnoise = fastnoise
        self._savemodel = savemodel
        self._parallel = parallel
        self._usemask = usemask

        ###
        ### 'self._mask' always contains the mask as supplied by the user while 'self._effective_mask' is
        ### the mask currently in play as interactive clean progresses. When the user has supplied a mask,
        ### it should be the same as 'self._mask' but when the mask is managed internally by iclean/gclean
        ### the two will diverge.
        ###
        self._mask = mask
        self.global_imdict = ImagingDict()
        self.current_imdict = ImagingDict()
        self._major_done = 0
        self.hasit = False # Convergence flag
        self.stopdescription = '' # Convergence flag
        self._convergence_result = (None,None,None,None,None,{ 'chan': None, 'major': None })
        #                           ^^^^ ^^^^ ^^^^ ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^----->>> convergence info
        #                              |    | |     |    +----->>> Number of global iterations remaining for current run (niterleft)
        #                              |    | |     +---------->>> Number of major cycles remaining for current run (nmajorleft)
        #                              |    | +---------------->>> major cycles done for current run (nmajordone)
        #                              |    +------------------>>> tclean stopcode
        #                              +----------------------->>> tclean stopdescription

        # Convert threshold from string to float, interpreting units.
        # XXX : We should ideally use quantities, but we are trying to
        # stick to "public API" funtions inside _gclean
        self._threshold_to_float()


    def __add_per_major_items( self, tclean_ret, major_ret, chan_ret ):
        '''Add meta-data about the whole major cycle, including 'cyclethreshold'
        '''

        if 'cyclethreshold' in tclean_ret:

            rdict = dict( major=dict( cyclethreshold=[tclean_ret['cyclethreshold']] if major_ret is None else (major_ret['cyclethreshold'] + [tclean_ret['cyclethreshold']]) ),
                         chan=chan_ret )
        else:
            rdict = dict( major=dict( cyclethreshold=major_ret['cyclethreshold'].append(tclean_ret['cyclethreshold']) ),
                         chan=chan_ret )

        return rdict


    def _calc_deconv_controls(self, imdict, niterleft=0, threshold=0, cycleniter=-1):
        """
        Calculate cycleniter and cyclethreshold for deconvolution.
        """

        use_cycleniter = niterleft  #niter - imdict.returndict['iterdone']

        if cycleniter > -1 : # User is forcing this number
            use_cycleniter = min(cycleniter, use_cycleniter)

        psffrac = imdict.returndict['maxpsfsidelobe'] * self._cyclefactor
        psffrac = max(psffrac, self._minpsffraction)
        psffrac = min(psffrac, self._maxpsffraction)

        # TODO : This assumes the default field (i.e., field=0);
        # This won't work for multiple fields.
        cyclethreshold = psffrac * imdict.get_peakres()
        cyclethreshold = max(cyclethreshold, threshold)

        return int(use_cycleniter), cyclethreshold


    def __update_convergence(self):
        """
        Accumulates the per-channel/stokes summaryminor keys across all major cycle calls so far.

        The "iterDone" key will be replaced with "iterations", and for the "iterations" key,
        the value in the returned cummulative record will be a rolling sum of iterations done
        for tclean calls so far, one value per minor cycle.
        For example, if there have been two clean calls, and in the first call channel 0 had
        [1] iteration in 1 minor cycle, and for the second call channel 0 had [6, 10, 9, 1]
        iterations in 4 minor cycles), then the resultant "iterations" key for channel 0 would be:
        [1, 7, 17, 26, 27]
        """

        keys = ['modelFlux', 'iterDone', 'peakRes', 'stopCode', 'cycleThresh']

        # Grab tuples of keys of interest
        outrec = {}
        for nn in range(self.global_imdict.nchan):
            outrec[nn] = {}
            for ss in range(self.global_imdict.nstokes):
                outrec[nn][ss] = {}
                for key in keys:
                    # Replace iterDone with iterations
                    if key == 'iterDone':
                        # Maintain cumulative sum of iterations per entry
                        outrec[nn][ss]['iterations'] = np.cumsum(self.global_imdict.get_key(key, stokes=ss, chan=nn))
                        # Replace iterDone with iterations
                        #outrec[nn][ss]['iterations'] =  self.global_imdict.get_key(key, stokes=ss, chan=nn)
                    else:
                        outrec[nn][ss][key] = self.global_imdict.get_key(key, stokes=ss, chan=nn)

        return outrec


    def __next__( self ):
        """ Runs tclean and returns the (stopcode, convergence result) when executed with the python builtin next() function.

        The returned convergence result is a nested dictionary:
        {
            channel id: {
                stokes id: {
                    summary key: [values, one per minor cycle]
                },
            },
        }

        See also: gclean.__update_convergence(...)
        """

        #                      vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv------------>>> is done to produce an initial dirty image
        if self._niter < 1 and self._convergence_result[2] is not None:
            self._convergence_result = ( f'nothing to run, niter == {self._niter}',
                                         self._convergence_result[1],
                                         self._major_done,
                                         self._nmajor,
                                         self._niter,
                                         self._convergence_result[3] )
            return self._convergence_result
        else:
            ### CALL SEQUENCE:
            ###       tclean(niter=0),deconvolve(niter=0),tclean(niter=100),deconvolve(niter=0),tclean(niter=100),tclean(niter=0,restoration=True)
            self._exe_cmds_per_iter.append(0)
            if self._convergence_result[1] is None:

                # initial call to tclean(...) creates the initial dirty image with niter=0
                tclean_ret = self._tclean( vis=self._vis, mask=self._mask, imagename=self._imagename, imsize=self._imsize, cell=self._cell,
                                           phasecenter=self._phasecenter, stokes=self._stokes, startmodel=self._startmodel, specmode=self._specmode,
                                           reffreq=self._reffreq, gridder=self._gridder, wprojplanes=self._wprojplanes, mosweight=self._mosweight,
                                           psterm=self._psterm, wbawp=self._wbawp, conjbeams=self._conjbeams, usepointing=self._usepointing,
                                           interpolation=self._interpolation, perchanweightdensity=self._perchanweightdensity,
                                           nchan=self._nchan, start=self._start, width=self._width, veltype=self._veltype, restfreq=self._restfreq,
                                           outframe=self._outframe, pointingoffsetsigdev=self._pointingoffsetsigdev, pblimit=self._pblimit,
                                           deconvolver=self._deconvolver, smallscalebias=self._smallscalebias, cyclefactor=self._cyclefactor,
                                           scales=self._scales, restoringbeam=self._restoringbeam, pbcor=self._pbcor, nterms=self._nterms,
                                           field=self._field, spw=self._spw, timerange=self._timerange, uvrange=self._uvrange, antenna=self._antenna,
                                           scan=self._scan, observation=self._observation, intent=self._intent, datacolumn=self._datacolumn,
                                           weighting=self._weighting, robust=self._robust, npixels=self._npixels, interactive=False, niter=0,
                                           gain=self._gain, calcres=True, calcpsf=True, restoration=False, parallel=self._parallel, fullsummary=True)


                deconv_ret = self._deconvolve(imagename=self._imagename, startmodel=self._startmodel,
                                              deconvolver=self._deconvolver, restoration=False,
                                              threshold=self._threshold, niter=0,
                                              nsigma=self._nsigma, fullsummary=True, fastnoise=self._fastnoise, usemask=self._usemask,
                                              mask=self._mask, noisethreshold=self._noisethreshold)

                self.current_imdict.returndict = self.current_imdict.merge(tclean_ret, deconv_ret)
                self.global_imdict.returndict = self.current_imdict.returndict

                # TODO : Add a standalone module to calculate the max PSF sidelobe.
                # Add it into deconv_ret at this point. The function can live inside imager_return_dict.py

                ## Initial call where niterleft and nmajorleft are same as original input values.
                self.hasit, self.stopdescription = self.global_imdict.has_converged(self._niter, self.current_imdict.get_key('threshold'), self._nmajor)

                self.current_imdict.returndict['stopcode'] = self.hasit
                self.current_imdict.returndict['stopDescription'] = self.stopdescription
                self._major_done = 0
            else:
                # Reset convergence every time, since we return control to the GUI after a single major cycle
                #self.hasit = 0
                #self.stopdescription = ''
                self.current_imdict.returndict['iterdone'] = 0.

                ### Check before doing the next round....
                self.hasit, self.stopdescription = self.global_imdict.has_converged(self._niter, self.global_imdict.get_key('threshold'), self._nmajor)

                #self.global_imdict.returndict['stopcode'] = self.hasit
                #self.global_imdict.returndict['stopDescription'] = self.stopdescription

                #self.current_imdict.returndict['stopcode'] = self.hasit
                #self.current_imdict.returndict['stopDescription'] = self.stopdescription

                if self.hasit ==0 :
                    use_cycleniter, cyclethreshold = self._calc_deconv_controls(self.current_imdict, self._niter, self._threshold, self._cycleniter)

                    # Run the minor cycle
                    deconv_ret = self._deconvolve(imagename=self._imagename, startmodel=self._startmodel,
                                              deconvolver=self._deconvolver, restoration=False,
                                              threshold=cyclethreshold, niter=use_cycleniter, gain=self._gain, usemask=self._usemask,
                                              nsigma=self._nsigma, fullsummary=True, fastnoise=self._fastnoise, noisethreshold=self._noisethreshold)

                    # Run the major cycle
                    tclean_ret = self._tclean( vis=self._vis, imagename=self._imagename, imsize=self._imsize, cell=self._cell,
                                           phasecenter=self._phasecenter, stokes=self._stokes, specmode=self._specmode, reffreq=self._reffreq,
                                           gridder=self._gridder, wprojplanes=self._wprojplanes, mosweight=self._mosweight, psterm=self._psterm,
                                           wbawp=self._wbawp, conjbeams=self._conjbeams, usepointing=self._usepointing, interpolation=self._interpolation,
                                           perchanweightdensity=self._perchanweightdensity, nchan=self._nchan, start=self._start,
                                           width=self._width, veltype=self._veltype, restfreq=self._restfreq, outframe=self._outframe,
                                           pointingoffsetsigdev=self._pointingoffsetsigdev, pblimit=self._pblimit, deconvolver=self._deconvolver,
                                           smallscalebias=self._smallscalebias, cyclefactor=self._cyclefactor, scales=self._scales,
                                           restoringbeam=self._restoringbeam, pbcor=self._pbcor, nterms=self._nterms, field=self._field,
                                           spw=self._spw, timerange=self._timerange, uvrange=self._uvrange, antenna=self._antenna,
                                           scan=self._scan, observation=self._observation, intent=self._intent, datacolumn=self._datacolumn,
                                           weighting=self._weighting, robust=self._robust, npixels=self._npixels, interactive=False,
                                           niter=0, restart=True, calcpsf=False, calcres=True, restoration=False, threshold=self._threshold,
                                           nsigma=self._nsigma, cycleniter=self._cycleniter, nmajor=1, gain=self._gain,
                                           sidelobethreshold=self._sidelobethreshold, noisethreshold=self._noisethreshold,
                                           lownoisethreshold=self._lownoisethreshold, negativethreshold=self._negativethreshold,
                                           minbeamfrac=self._minbeamfrac, growiterations=self._growiterations, dogrowprune=self._dogrowprune,
                                           minpercentchange=self._minpercentchange, fastnoise=self._fastnoise, savemodel=self._savemodel, maxpsffraction=self._maxpsffraction,
                                           minpsffraction=self._minpsffraction, parallel=self._parallel, fullsummary=True )

#                    print("DECONV RET")
#                    print(deconv_ret)
#                    print("TCLEAN RET")
#                    print(tclean_ret)

                    # Replace return dict with new return dict
                    # The order of the dicts into merge is important.
                    self.current_imdict.returndict = self.current_imdict.merge(tclean_ret, deconv_ret)
                    # Append new return dict to global return dict
                    self.global_imdict.returndict = self.global_imdict.concat(self.global_imdict.returndict, self.current_imdict.returndict)
                    self._major_done = self.current_imdict.returndict['nmajordone']

                    ## Decrement count for the major cycle just done...
                    self.__decrement_counts()

                    # Use global imdict for convergence check
                    self.hasit, self.stopdescription = self.global_imdict.has_converged(self._niter, self.global_imdict.get_key('threshold'), self._nmajor)


                self.global_imdict.returndict['stopcode'] = self.hasit
                self.global_imdict.returndict['stopDescription'] = self.stopdescription

                if not self.hasit:
                    # If we haven't converged, run deconvolve to update the mask
                    self._deconvolve(imagename=self._imagename, niter=0, deconvolver=self._deconvolver, usemask=self._usemask, restoration=False)

            if len(self.global_imdict.returndict) > 0 and 'summaryminor' in self.global_imdict.returndict and sum(map(len,self.global_imdict.returndict['summaryminor'].values())) > 0:
                # self.current_imdict only contains the latest tclean/deconvolve results
                # Passing in self.global_imdict will pull out the cumulative results everytime, breaking the convergence plot.
                self._convergence_result = ( self.global_imdict.returndict['stopDescription'] if 'stopDescription' in self.global_imdict.returndict else '',
                                             self.global_imdict.returndict['stopcode'] if 'stopcode' in self.global_imdict.returndict else 0,
                                             self._major_done,
                                             self._nmajor,
                                             self._niter,
                                             self.__add_per_major_items( self.global_imdict.returndict,
                                                                         self._convergence_result[5]['major'],
                                                                         self.__update_convergence()))
            else:
                self._convergence_result = ( f'tclean returned an empty result',
                                             self._convergence_result[1],
                                             self._major_done,
                                             self._nmajor,
                                             self._niter,
                                             self._convergence_result[5] )

            return self._convergence_result

    def __decrement_counts( self ):
        ## Update niterleft and nmajorleft now.
        if self.hasit == 0:  ##If not yet converged.
            if self._nmajor != -1:   ## If -1, don't touch it.
                self._nmajor = self._nmajor - 1
            self._niter = self._niter - self.current_imdict.get_key('iterdone')
            if self._niter<0:  ## This can happen when we're counting niter across channels in a single minor cycle set, and it crosses the total. 
                self._niter=0
        else:
            return  ##If convergence has been reached, don't try to decrement further.


    def __reflect_stop( self ):
        ## if python wasn't hacky, you would be able to try/except/raise in lambda
        time.sleep(1)
        try:
            return self.__next__( )
        except StopIteration:
            raise StopAsyncIteration

    async def __anext__( self ):
        ### asyncio.run cannot be used here because this is called
        ### within an asyncio loop...
        loop = asyncio.get_event_loop( )
        result = await loop.run_in_executor( None, self.__reflect_stop )
        return result

    def __iter__( self ):
        return self

    def __aiter__( self ):
        return self

    def __split_filename( self, path ):
        return os.path.splitext(os.path.basename(path))

    def __default_mask_name( self ):
        imgparts = self.__split_filename( self._imagename )
        return f'''{imgparts[0]}.mask'''

    def mask(self):
        return self.__default_mask_name() if self._mask == '' else self._mask

    def reset(self):
        #if not self._finalized:
        #    raise RuntimeError('attempt to reset a gclean run that has not been finalized')
        self._finalized = False
        self._convergence_result = ( None,
                                     self._convergence_result[1],
                                     self._major_done,
                                     self._nmajor,
                                     self._niter,
                                     self._convergence_result[5] )

    def restore(self):
        """ Restores the final image, and returns a path to the restored image. """
        deconv_ret = self._deconvolve(imagename=self._imagename, startmodel=self._startmodel,
                                      deconvolver=self._deconvolver, restoration=True,
                                      threshold=self._threshold, niter=0, gain=self._gain,
                                      nsigma=self._nsigma, fullsummary=True, fastnoise=self._fastnoise, usemask=self._usemask,
                                      noisethreshold=self._noisethreshold)

        return { "image": f"{self._imagename}.image" }

    def has_next(self):
        return not self._finalized
