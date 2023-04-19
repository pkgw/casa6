##########################################################################
#
# Run the tests as described in
# # https://open-confluence.nrao.edu/pages/viewpage.action?spaceKey=CASA&title=Requirements+for+VLASS+Imaging+Pipeline+Stakeholders+Tests
# There wasn't a great definition of the work to be done, so much of it is being interpretted by me (BGB).
#
##########################################################################
#
# The basic idea:
#  Compare values from the most recent build of casa to other known values. Check for differences.
#
# Values to compare current versions against:
#  1. on-axis      ---   "ground truth", I think from the (fluxscale?) task
#  2. CASA 6.2.1   ---   pipeline approved version of casa
#
# What's good enough?
# Empirical tolerance:
#  flux density: 5% goal, 10% ok
#  spectral index: 0.1 goal, 0.2 ok
#  reference: https://drive.google.com/file/d/1zw6UeDEoXoxM05oFg3rir0hrCMEJMxkH/view and https://open-confluence.nrao.edu/display/VLASS/Updated+VLASS+survey+science+requirements+and+parameters
# Noise floor adjusted tolerance:
#  flux density & F_nu: (2*rms) / max(|expected|,2*rms)
#  spectral index: truth*sqrt( (2*rms0/truth0)^2 + (2*rms1/truth1)^2 )
#  reference: https://casadocs.readthedocs.io/en/latest/notebooks/synthesis_imaging.html#Options-in-CASA-for-wideband-imaging --> Search for "Calculating Error in Spectral Index"
#  reference: Eqn 39 of https://www.aanda.org/index.php?option=com_article&access=doi&doi=10.1051/0004-6361/201117104&Itemid=129#S29
#
# Images:
#  J1302
#  J1927
#
# Image processing methods:
#  mtmfs         ---   "stokes I"
#  awproject     ---   "stokes I"
#  mosaic cube   ---   "cube"
#  mosaic QL     ---   "QL"
# Notes:
#  mtmfs methods referred to as "mosaic"
#  "mosaic cube" is actual a vlass "coarse cube" (manual mtmfs cube via multiple tclean runs for select spws)
#  "ql" is short for "quick look", a faster way of processing images for non-scientific results for vlass
#
# Values to be verified:                                       Compare against:
#  Stokes I (and QL?):
#   a. tt0:                                                    6.2.1, on-axis
#   b. tt1:                                                    6.2.1, on-axis
#   c. alpha:                                                  6.2.1, on-axis
#   d. beamsize comparison:                                    6.2.1
#  Stokes I:
#   e. Confirm presence of model column in resultant MS
#  Stokes I and Cube:
#   f. Runtimes not significantly different relative to previous runs
#  Cube:
#   g. Fit F_nu0 and Alpha from three cube planes and compare: 6.2.1, on-axis
#   h. IQUV flux densities of all three spws:                  6.2.1
#   i. IQUV flux densities of all three spws:                  on-axis measurements
#   j. Beam of all three spws:                                 6.2.1
#  QL:
#   k. flux density of Calibrator source:                      6.2.1, on-axis
#  All:
#   l. Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
# Note on (i): on-axis values derived from images of calibrator using standard gridder of pointed data on calibrator
#              ^-- does this mean these are also "ground truth" values? BGB211222
#
# Other goals:
#  AUX1: run each test as its own test script, eg "VLASS_mosaic_stakeholder_test_script.py"
#   not done (is this actually necessary? what is the goal of these scripts such that individual files are needed?)
#
##########################################################################
#
#  Tests
#About Two Data sets (based on John Tobin's comment in CAS-13832)
# J1302: southern source (-10.51.16.73), has a flat spectrum with a little polarization
# J1927: northern source (+61.17.32.898), has a steeper spectrum and polarized emission

#J1302 Tests
#1. mtmfs: Should match values for "Stokes I" in the "Values to be compared"
#vis:'J1302_12field.ms' gridder:'mosaic'
#testname: test_j1302_mosaic
#
#2. awproject: Should match values for "Stokes I" in the "Values to be compared"
#vis:'J1302_12field.ms' gridder:'awproject'
#testname: test_j1302_awproject
#
#3. mosaic cube: Should match values for "Cube" in the "Values to be compared"
#vis:'J1302_12field_cubedata.ms' gridder:'mosaic'
#testname: test_j1302_mosaic_cube
#
#4. QL: Should match values for "QL" in the "Values to be compared"
#vis:'J1302_12field.ms' gridder:'ql'
#testname: test_j1302_ql
#
#
#
#J1927 Tests
#5. mtmfs: Should match values for "Stokes I" in the "Values to be compared"
#vis:'J1927_12field.ms', gridder:'mosaic'
#testname: test_j1927_mosaic
#
#6. mosaic cube: Should match values for "Cube" in the "Values to be compared"
#vis:'J1927_12field_cubedata.ms', gridder:'mosaic'
#testname: test_j1927_mosaic_cube
#
#7. QL: Should match values for "QL" in the "Values to be compared"
#vis:'J1927_12field.ms', gridder:'ql'
#testname: test_j1927_ql
#
##########################################################################

import os
import unittest
import numpy as np
import shutil
import glob
from datetime import datetime

from casatools import imager
from casatasks import casalog, impbcor, imdev, imhead, imsubimage, imstat, immath
from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
from casatestutils.stakeholder import StkUnitTest
from casatestutils import generate_weblog
from casatestutils import add_to_dict
from casatestutils import stats_dict

quick_test = False if ('QUICK_TEST' not in os.environ) else os.environ['QUICK_TEST']
quick_test = True if str(quick_test).lower() in ['1', 'true'] else False
use_partial_results = False if ('USE_PARTIAL_RESULTS' not in os.environ) else os.environ['USE_PARTIAL_RESULTS']
use_partial_results = True if str(use_partial_results).lower() in ['1', 'true'] else False
cache_partial_results = False if ('CACHE_PARTIAL_RESULTS' not in os.environ) else os.environ['CACHE_PARTIAL_RESULTS']
cache_partial_results = True if str(cache_partial_results).lower() in ['1', 'true'] else False
if (quick_test):
    casalog.post("QUICK_TEST env variable found\nRunning tests with reduced image sizes", "INFO")
else:
    casalog.post("QUICK_TEST env variable not found or is false\nRunning tests with full image sizes", "INFO")
casalog.post(f"USE_PARTIAL_RESULTS: {use_partial_results}", "INFO")


test_dict = {}

##############################################
##############################################
class test_j1302(StkUnitTest):

    def setUp(self):
        super().setUp()
        self.vis = 'J1302-12fields.ms'
        self.phasecenter = '13:03:13.874 -10.51.16.73'
        self.im = imager()
        self.parallel = False
        if ParallelTaskHelper.isMPIEnabled():
            self.parallel = True


    def tearDown(self):
        generate_weblog("tclean_VLASS_1v2_pipeline", test_dict)
        super().tearDown()
        if not cache_partial_results:
            self.delData() 
        
    # move delData in StkUnitTest here
    def delData(self):
        """ Clean up generated data. """
        # clean up files listed in self.vis, self.imgs, and
        # self.teardown_files
        del_files = []
        if self.vis != "" and os.path.exists(self.vis):
            del_files.append(self.vis)
        for img in self.imgs:
            img_files = glob.glob(img+'*')
            del_files += img_files
        for teardown_file in self.teardown_files:
            if teardown_file in del_files:
                continue
            if not os.path.exists(teardown_file):
                continue
            del_files.append(teardown_file)

        # don't delete weblogs at the end (done in setUpClass instead)
        keep_files = list(filter(lambda f: f.endswith(".png") or f.endswith(".html"), del_files))
        del_files  = list(filter(lambda f: f not in keep_files, del_files))

        # delete the del_files
        for f in del_files:
            self.del_file_or_dir(f)


    # Test 1
    @stats_dict(test_dict)
    def test_j1302_mtmfs(self):
        """ [j1302] test_j1302_mtmfs """
        ######################################################################################
        # Should match values for "Stokes I" in the "Values to be compared"
        ######################################################################################
        # "noncube" to allow give this test a unique prefix, for running with runtest
        # For example: runtest.py -v test_vlass_1v2.py[test_j1302_mtmfs]

        ##############################################
        # %% Set local vars [test_j1302_mtmfs] start @
        ##############################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        test_name = self._testMethodName
        data_path_dir = 'J1302/Stakeholder-test-mosaic-data'
        #img0 = 'J1302_iter2'
        basetname = 'J1302_mtmfs'
        img0 = f"{basetname}_iter2"
        masks = ['secondmask.mask', 'QLcatmask.mask']
        quick_masks = ['secondmask_1000.mask', 'QLcatmask_1000.mask']

        if not quick_test:
            imsize=4000
            self.prepData(self.vis, data_path_dir, *masks, partial_results_dirname="partial_results_test_j1302_mtmfs")
            for i in range(len(masks)):
                os.system(f"mv {masks[i]} {basetname}_{masks[i]}")
        else:
            imsize=1000
            self.prepData(self.vis, data_path_dir, *quick_masks, partial_results_dirname="partial_results_test_j1302_mtmfs")
            for i in range(len(masks)):
                os.system(f"mv {quick_masks[i]} {basetname}_{masks[i]}")
        self.teardown_files += [f'{basetname}_{x}' for x in masks]

        #spw = ''
        rms = [0.00017975829898762892, 0.0013099727978948515] # tt0, tt1 noise floor as measured from a full-scale image run, Range: [700,800],[3300,1900]
        starttime = datetime.now()

        #############################################
        # %% Set local vars [test_j1302_mtmfs] end  @
        # %% Prepare masks [test_j1302_mtmfs] start @
        #############################################

        # combine first and 2nd order masks
        if not use_partial_results:
            immath(imagename=[f'{basetname}_secondmask.mask',f'{basetname}_QLcatmask.mask'],expr='IM0+IM1',outfile=f'{basetname}_sum_of_masks.mask')
            self.im.mask(image=f'{basetname}_sum_of_masks.mask',mask=f'{basetname}_combined.mask',threshold=0.5)
            self.teardown_files += [f'{basetname}_sum_of_masks.mask', f'{basetname}_combined.mask']

        ###########################################
        # %% Prepare masks [test_j1302_mtmfs] end @
        # %% Run tclean [test_j1302_mtmfs] start  @
        ###########################################

        records = []
        tstobj = self
        common_args = {'vis':'J1302-12fields.ms', 'uvrange':'<12km', 'intent':'OBSERVE_TARGET#UNSPECIFIED', \
                       'imagename':img0, 'imsize':imsize, 'cell':'0.6arcsec', 'phasecenter':'13:03:13.874 -10.51.16.73', \
                       'spw':'', 'reffreq':'3.0GHz', \
                       'gridder':'mosaic', 'conjbeams':False, 'mosweight':False, 'pblimit':0.1, \
                       'deconvolver':'mtmfs', 'nterms':2, 'smallscalebias':0.4, \
                       'weighting':'briggs', 'robust':1.0,\
                       'threshold':0.0, 'cyclefactor':3.0, 'restart':True, 'interactive':False, 'parallel':self.parallel}

        def run_tclean(niter=None, datacolumn='corrected', mask='', usemask='user', pbmask=0.0, nsigma=2.0, cycleniter=500, scales=[0],
              calcres=False, calcpsf=False, savemodel='none', common_args=common_args):
            params = locals()
            finalparams = self.filter_runtclean_parameters(params, common_args)
            # add images to self.imgs for clean-up
            if ('imagename' in finalparams):
                img = finalparams['imagename']
                if (img not in self.imgs):
                    self.imgs.append(img)
            return tclean(**finalparams)

            #params = {k:params[k] for k in filter(lambda x: x not in ['tstobj', 'records'], params.keys())}
            #records.append( tstobj.run_tclean(**params) )
            #return records[-1]

        # These values are lifted from logs with John Tobin's original test scripts on the CAS-12427 ticket.
        # They are only used as a comparison, to make sure we use all the same exact parameter vlaues.
        #script_pars_vals_0 = self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False)

        #script_pars_vals_1 = self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)
        #script_pars_vals_2 = self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='J1302_iter2', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='modelcolumn', calcres=False, calcpsf=False, parallel=False)
        #script_pars_vals_3 = self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)
        #script_pars_vals_4 = self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)

        ##### main test execution begin #####
        # initialize iter2, no cleaning
        #run_tclean( niter=0,     datacolumn='corrected', calcres=True, calcpsf=True,                         compare_tclean_pars=script_pars_vals_0 )
        run_tclean( niter=0, calcres=True, calcpsf=True)

        # resume iter2 with QL mask
        #run_tclean( niter=20000, datacolumn='corrected', mask=f"{basetname}_QLcatmask.mask", nsigma=3.0, scales=[0,5,12], compare_tclean_pars=script_pars_vals_1 )
        run_tclean( niter=20000, mask=f"{basetname}_QLcatmask.mask", nsigma=3.0, scales=[0,5,12] ) 

        # save model column, doesn't happen here in acutal VLASS pipeline, but makes sure functionality works.
        #run_tclean( niter=0,     datacolumn='data',      savemodel='modelcolumn',                            compare_tclean_pars=script_pars_vals_2 )
        run_tclean( niter=0,     datacolumn='data',      savemodel='modelcolumn' )

        # resume iter2 with combined mask, remove old mask first, pass new mask as parameter
        os.system(f"rm -rf {img0}.mask")
        #run_tclean( niter=20000, datacolumn='corrected', mask=f"{basetname}_combined.mask",  nsigma=3.0, scales=[0,5,12], compare_tclean_pars=script_pars_vals_3 )
        run_tclean( niter=20000,  mask=f"{basetname}_combined.mask",  nsigma=3.0, scales=[0,5,12] )

        # resume iter2 with pbmask, removed old mask first then specify pbmask in resumption of tclean
        os.system(f"rm -rf {img0}.mask")
        #run_tclean( niter=20000, datacolumn='corrected', usemask='pb', mask="", pbmask=0.4,   nsigma=4.5, scales=[0,5,12], cycleniter=100, compare_tclean_pars=script_pars_vals_4 )
        run_tclean( niter=20000, usemask='pb', mask="", pbmask=0.4,   nsigma=4.5, scales=[0,5,12], cycleniter=100 )
        ##### main test execution end #####

        #######################################################
        # %% Run tclean [test_j1302_mtmfs] end                @
        # %% Compare Expected Values [test_j1302_mtmfs] start @
        #######################################################

        stats621 = {
            'full': {
                'tt0': 0.31962126,
                'tt1': 0.01928869,
                'alpha': 0.06034857,
                'beam': { 'maj': 3.15002894, 'min': 2.59177852, 'pos': 11.41288376 },
                #'runtime': 2413
            },
            '1/4': {
                'tt0': 0.29668888,
                'tt1': -0.00507728,
                'alpha': -0.01711315,
                'beam': { 'maj': 3.13147807, 'min': 2.56502342, 'pos': 10.92203331 },
                #'runtime': 417
            }
        }
        stats621 = stats621['full'] if (not quick_test) else stats621[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        # N/A for this test

        halfsize   = round(imsize / 4000 * 2000)
        box        = f"{halfsize},{halfsize},{halfsize},{halfsize}"
        tt0stats   = imstat(imagename=img0+'.image.tt0', box=box)
        tt1stats   = imstat(imagename=img0+'.image.tt1', box=box)
        alphastats = imstat(imagename=img0+'.alpha', box=box)
        curr_stats    = np.squeeze(np.array([ tt0stats['max'], tt1stats['max'], alphastats['max']]))
        onaxis_stats  = np.array([            0.3337,         -0.01588,        -0.0476])
        casa621_stats = np.array([            stats621['tt0'], stats621['tt1'], stats621['alpha']])

        # (a) tt0 vs 6.2.1, on-axis
        success0, report0 = self.check_metrics_flux(curr_stats[0], onaxis_stats[0],  valname="Frac Diff tt0 vs. on-axis", rms_or_std=rms[0])
        success1, report1 = self.check_metrics_flux(curr_stats[0], casa621_stats[0], valname="Frac Diff tt0 vs. 6.2.1 image", rms_or_std=rms[0])

        # (b) tt1 vs 6.2.1, on-axis
        success2, report2 = self.check_metrics_flux(curr_stats[1], onaxis_stats[1],  valname="Frac Diff tt1 vs. on-axis", rms_or_std=rms[1])
        success3, report3 = self.check_metrics_flux(curr_stats[1], casa621_stats[1], valname="Frac Diff tt1 vs. 6.2.1 image", rms_or_std=rms[1])

        # (c) alpha images
        success4, report4 = self.check_metrics_alpha(curr_stats[2], onaxis_stats[2],  onaxis_stats[0], onaxis_stats[1], valname="Abs Diff alpha vs. on-axis", rmss_or_stds=rms)
        success5, report5 = self.check_metrics_alpha(curr_stats[2], casa621_stats[2], casa621_stats[0], casa621_stats[1], valname="Abs Diff alpha vs. 6.2.1 image", rmss_or_stds=rms)

        # (d) beamsize comparison vs 6.2.1
        restbeam          = imhead(img0+'.image.tt0')['restoringbeam']
        beamstats_curr    = np.array([restbeam['major']['value'], restbeam['minor']['value'], restbeam['positionangle']['value']])
        beamstats_621     = np.array([stats621['beam']['maj'],    stats621['beam']['min'],    stats621['beam']['pos']])
        success6, report6 = self.check_fracdiff(beamstats_curr, beamstats_621, valname="Frac Diff Maj, Min, PA vs 6.2.1")

        # (e) Confirm presence of model column in resultant MS
        success7, report7 = self.check_column_exists("MODEL_DATA")

        #####################################################
        # %% Compare Expected Values [test_j1302_mtmfs] end @
        # %% Generate Images [test_j1302_mtmfs] start       @
        #####################################################

        #self.mom8_creator(image=img0+'.image.tt0', range_list=[0.003, 0.32], imgname="j1302_mtmfs_tt0")
        #self.mom8_creator(image=img0+'.image.tt1', range_list=[0.001, 0.04], imgname="j1302_mtmfs_tt1")
        self.mom8_creator(image=img0+'.image.tt0', range_list=[-0.3, 0.32], imgname="j1302_mtmfs_tt0")
        self.mom8_creator(image=img0+'.image.tt1', range_list=[-0.1, 0.04], imgname="j1302_mtmfs_tt1")

        #############################################
        # %% Generate Images [test_j1302_mtmfs] end @
        #############################################

        # save results for future analysis
        np.save(self.id()+'.tt0tt1alpha.npy', curr_stats)
        np.save(self.id()+'.beamstats.npy', beamstats_curr)
        np.save(self.id()+'.tcleanrecs.npy', records)

        # (f) Runtimes not significantly different relative to previous runs
        # if ParallelTaskHelper.isMPIEnabled():
        #     # runtime with MPI -n 8
        #     successr, reportr = self.check_runtime(starttime, 5373)
        # else:
        #     successr, reportr = self.check_runtime(starttime, stats621['runtime'])

        # TODO start obeying on-axis once on-axis passes
        report  = "".join([report1, report3, report4, report5, report6, report7])
        success = success1 and success3 and success4 and success5 and success6 and success7 and self.th.check_final(report)
        # report  = "".join([report0, report1, report2, report3, report4, report5, report6, report7])
        # success = success0 and success1 and success2 and success3 and success4 and success5 and success6 and success7 and self.th.check_final(report)

        add_to_dict(self, output = test_dict, dataset = "J1302-12fields.ms")
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = self.mom8_images
        test_dict[test_name]['taskcall'] = self.clean_taskcall(test_dict[test_name]['taskcall'], locals())

        self.assertTrue(success, msg=self.th.extract_failing_lines(report))

    # Test 2
    @stats_dict(test_dict)
    def test_j1302_awproject(self):
        """ [j1302] test_j1302_awproject """
        ######################################################################################
        # Should match values for "Stokes I" in the "Values to be compared"
        ######################################################################################

        ##################################################
        # %% Set local vars [test_j1302_awproject] start @
        ##################################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        test_name = self._testMethodName
        data_path_dir  = 'J1302/Stakeholder-test-awproject-data'
        cf0 = 'J1302_iter0d'
        cf1 = 'J1302_iter2'
        #img0 = 'J1302_iter0d'
        #img1 = 'J1302_iter2'
        basetname = 'J1302_awproject'
        img0 = f'{basetname}_iter0d'
        img1 = f'{basetname}_iter2'
        cache0name, cache1name = "cache0d.cf", "cache2.cf"
        masks = ['secondmask.mask', 'QLcatmask.mask']
        quick_masks = ['secondmask_1312.mask', 'QLcatmask_1312.mask']

        if not quick_test:
            imsize=5250
            #self.prepData(self.vis, data_path_dir, f"cfcache/{img0}.cf", f"cfcache/{img1}.cf", *masks, partial_results_dirname="partial_results_test_j1302_awproject")
            self.prepData(self.vis, data_path_dir, f"cfcache/{cf0}.cf", f"cfcache/{cf1}.cf", *masks, partial_results_dirname="partial_results_test_j1302_awproject")
            #os.system(f"mv cfcache/{img0}.cf {cache0name}")
            #os.system(f"mv cfcache/{img1}.cf {cache1name}")
            os.system(f"mv cfcache/{cf0}.cf {cache0name}")
            os.system(f"mv cfcache/{cf1}.cf {cache1name}")
            for i in range(len(masks)):
                os.system(f"mv {masks[i]} {basetname}_{masks[i]}")
        else:
            imsize=1312
            #self.prepData(self.vis, data_path_dir, f"cfcache_quick1/{img0}.cf", f"cfcache_quick1/{img1}.cf", *quick_masks, partial_results_dirname="partial_results_test_j1302_awproject")
            #os.system(f"mv cfcache_quick1/{img0}.cf {cache0name}")
            #os.system(f"mv cfcache_quick1/{img1}.cf {cache1name}")
            self.prepData(self.vis, data_path_dir, f"cfcache_quick1/{cf0}.cf", f"cfcache_quick1/{cf1}.cf", *quick_masks, partial_results_dirname="partial_results_test_j1302_awproject")
            os.system(f"mv cfcache_quick1/{cf0}.cf {cache0name}")
            os.system(f"mv cfcache_quick1/{cf1}.cf {cache1name}")
            for i in range(len(masks)):
                os.system(f"mv {quick_masks[i]} {basetname}_{masks[i]}")
        self.teardown_files += [f'{basetname}_{x}' for x in masks]
        self.teardown_files += [cache0name, cache1name]

        if not quick_test:
            rms = [0.00025982361923319354, 0.00211483438886223] # tt0, tt1 noise floor as measured from a full-scale image run, Range: [1500,500],[3500,2000]
        else:
            rms = [0.00181162, 0.02026735] # tt0, tt1 noise floor as measured in the QUICK_TEST images with an annulus r=15 ~ 120 arcsec centered at the phase center
            # rms for the entire QUICK_TEST images
            #rms = [0.00143267, 0.01438686]
            
          
             
        starttime = datetime.now()

        ################################################
        # %% Set local vars [test_j1302_awproject] end @
        # %% Prepare masks [test_j1302_mtmfs] start    @
        ################################################

        # combine first and 2nd order masks
        if not use_partial_results:
            immath(imagename=[f'{basetname}_secondmask.mask',f'{basetname}_QLcatmask.mask'],expr='IM0+IM1',outfile=f'{basetname}_sum_of_masks.mask')
            self.im.mask(image=f'{basetname}_sum_of_masks.mask',mask=f'{basetname}_combined.mask',threshold=0.5)
            self.teardown_files += [f'{basetname}_sum_of_masks.mask', f'{basetname}_combined.mask', 
                                    f'{basetname}_secondmask.mask',f'{basetname}_QLcatmask.mask']

        ##############################################
        # %% Prepare masks [test_j1302_mtmfs] end    @
        # %% Run tclean [test_j1302_awproject] start @
        ##############################################

        records = []
        tstobj = self
        common_args = {'vis':'J1302-12fields.ms', 'uvrange':'<12km', 'intent':'OBSERVE_TARGET#UNSPECIFIED', \
                       'imsize':imsize, 'cell':'0.6arcsec', 'phasecenter':'13:03:13.874 -10.51.16.73', \
                       'spw':'', 'reffreq':'3.0GHz', \
                       'gridder':'awproject', 'wprojplanes':32, 'conjbeams':True, 'mosweight':False, 
                       'usepointing':True, 'rotatepastep':5.0, 'pointingoffsetsigdev':[300, 30], \
                       'pblimit':0.02, \
                       'deconvolver':'mtmfs', 'nterms':2, 'smallscalebias':0.4, \
                       'weighting':'briggs', 'robust':1.0,\
                       'threshold':0.0, 'cyclefactor':3.0, 'restart':True, 'interactive':False, 'parallel':self.parallel}

        def run_tclean(imagename=None, niter=None, datacolumn='corrected', mask='', usemask='user', pbmask=0.0, cfcache='', wbawp=True,  nsigma=2.0, cycleniter=5000, scales=[0],
              calcres=False, calcpsf=False, savemodel='none', common_args=common_args):
            params = locals()
            finalparams = self.filter_runtclean_parameters(params, common_args)
            # add images to self.imgs for clean-up
            if ('imagename' in finalparams):
                img = finalparams['imagename']
                if (img not in self.imgs):
                    self.imgs.append(img)
            return tclean(**finalparams)

        def replace_psf(old, new):
            """ Replaces [old] PSF image with [new] image. Clears parallel working directories."""
            for this_tt in ['tt0', 'tt1', 'tt2']:
                shutil.rmtree(old+'.psf.'+this_tt)
                shutil.copytree(new+'.psf.'+this_tt, old+'.psf.'+this_tt)

        # These values are lifted from logs with John Tobin's original test scripts on the CAS-12427 ticket.
        # They are only used as a comparison, to make sure we use all the same exact parameter vlaues.
        #script_pars_vals_0 = self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter0d', imsize=5250, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='awproject', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=32, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=False, conjbeams=True, cfcache='', usepointing=True, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[300, 30], pblimit=0.02, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=5000, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=True, parallel=False )
        #script_pars_vals_1 = self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=5250, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='awproject', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=32, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=True, cfcache='', usepointing=True, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[300, 30], pblimit=0.02, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=5000, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False )
        #script_pars_vals_2 = self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=5250, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='awproject', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=32, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=True, cfcache='', usepointing=True, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[300, 30], pblimit=0.02, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=3000, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False )
        #script_pars_vals_3 = self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='J1302_iter2', imsize=5250, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='awproject', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=32, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=True, cfcache='', usepointing=True, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[300, 30], pblimit=0.02, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=5000, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='modelcolumn', calcres=False, calcpsf=False, parallel=False )
        #script_pars_vals_4 = self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=5250, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='awproject', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=32, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=True, cfcache='', usepointing=True, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[300, 30], pblimit=0.02, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=3000, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False )
        #script_pars_vals_5 = self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=5250, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='awproject', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=32, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=True, cfcache='', usepointing=True, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[300, 30], pblimit=0.02, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False )

        # create robust psfs with wbawp=False using corrected column
        run_tclean(imagename=img0, cfcache=cache0name, niter=0, calcres=False, calcpsf=True, wbawp=False)

        # initialize iter2, no cleaning
        run_tclean(imagename=img1, cfcache=cache1name, niter=0, calcres=True, calcpsf=True )

        #  replace iter2 psf with no_WBAP
        replace_psf(img1, img0)

        #  resume iter2 with QL mask
        run_tclean(imagename=img1, cfcache=cache1name, niter=20000, scales=[0, 5, 12], nsigma=3.0, cycleniter=3000,
                   mask=f"{basetname}_QLcatmask.mask")

        # save model column, doesn't happen here in acutal VLASS pipeline, but makes sure functionality works.
        run_tclean(imagename=img1, cfcache=cache1name, niter=0, datacolumn='data', savemodel='modelcolumn')

        # resume iter2 with combined mask
        os.system(f"rm -rf {img1}.mask")
        run_tclean(imagename=img1, cfcache=cache1name, niter=20000, scales=[0, 5, 12], nsigma=3.0, cycleniter=3000,
                   mask=f"{basetname}_combined.mask")

        # resume iter2 with pbmask, removed old mask first then specify pbmask in resumption of tclean
        os.system(f"rm -rf {img1}.mask")
        run_tclean(imagename=img1, cfcache=cache1name, niter=20000, scales=[0, 5, 12], nsigma=4.5, cycleniter=500,
                   mask="", usemask='pb', pbmask=0.4)

        ###########################################################
        # %% Run tclean [test_j1302_awproject] end                @
        # %% Compare Expected Values [test_j1302_awproject] start @
        ###########################################################

        stats621 = {
            'full': {
                'tt0': 0.317386,
                'tt1': -0.01537329,
                'alpha': -0.04843721,
                'beam': { 'maj': 3.05415082, 'min': 2.49866414, 'pos': 10.82348919 },
                #'runtime': 21179
            },
            '1/4': {
                'tt0': 0.21766607,
                'tt1': -0.05709893,
                'alpha': -0.2623235,
                'beam': { 'maj': 3.07264543, 'min': 2.5105536,  'pos': 10.96454144 },
                #'runtime': 1993 # this is with cfcache
            }
        }
        stats621 = stats621['full'] if (not quick_test) else stats621[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        # N/A: no pbcore, rms, or subim images are created for this test

        halfsize   = round(imsize / 5250 * 2625)
        box        = f"{halfsize},{halfsize},{halfsize},{halfsize}"
        tt0stats=imstat(  imagename=img1+'.image.tt0',box=box)
        tt1stats=imstat(  imagename=img1+'.image.tt1',box=box)
        alphastats=imstat(imagename=img1+'.alpha',    box=box)
        curr_stats=np.squeeze(np.array([tt0stats['max'], tt1stats['max'], alphastats['max']]))
        onaxis_stats=np.array([         0.3337,         -0.01588,        -0.0476])
        casa621_stats=np.array([        stats621['tt0'], stats621['tt1'], stats621['alpha']])

        # (a) tt0 vs 6.2.1, on-axis
        success0, report0 = self.check_metrics_flux(curr_stats[0], onaxis_stats[0],  valname="Frac Diff tt0 vs. on-axis", rms_or_std=rms[0])
        success1, report1 = self.check_metrics_flux(curr_stats[0], casa621_stats[0], valname="Frac Diff tt0 vs. 6.2.1 image", rms_or_std=rms[0])

        # (b) tt1 vs 6.2.1, on-axis
        success2, report2 = self.check_metrics_flux(curr_stats[1], onaxis_stats[1],  valname="Frac Diff tt1 vs. on-axis", rms_or_std=rms[1])
        success3, report3 = self.check_metrics_flux(curr_stats[1], casa621_stats[1], valname="Frac Diff tt1 vs. 6.2.1 image", rms_or_std=rms[1])

        # (c) alpha images
        success4, report4 = self.check_metrics_alpha(curr_stats[2], onaxis_stats[2],  onaxis_stats[0], onaxis_stats[1], valname="Abs Diff alpha vs. on-axis", rmss_or_stds=rms)
        success5, report5 = self.check_metrics_alpha(curr_stats[2], casa621_stats[2], casa621_stats[0], casa621_stats[1], valname="Abs Diff alpha vs. 6.2.1 image", rmss_or_stds=rms)

        # (d) beamsize comparison vs 6.2.1
        restbeam          = imhead(img1+'.image.tt0')['restoringbeam']
        beamstats_curr    = np.array([restbeam['major']['value'], restbeam['minor']['value'], restbeam['positionangle']['value']])
        beamstats_621     = np.array([stats621['beam']['maj'],    stats621['beam']['min'],    stats621['beam']['pos']])
        success6, report6 = self.check_fracdiff(beamstats_curr, beamstats_621, valname="Frac Diff Maj, Min, PA vs 6.2.1")

        # (e) Confirm presence of model column in resultant MS
        success7, report7 = self.check_column_exists("MODEL_DATA")

        #########################################################
        # %% Compare Expected Values [test_j1302_awproject] end @
        # %% Generate Images [test_j1302_awproject] start       @
        #########################################################

        #self.mom8_creator(image=img1+'.image.tt0', scaling=-2, range_list=[0.002, 0.32], imgname="j1302_awproject_iter2_tt0")
        #self.mom8_creator(image=img1+'.image.tt1', scaling=-2, range_list=[0.007, 0.16], imgname="j1302_awproject_iter2_tt1")
        self.mom8_creator(image=img1+'.image.tt0', range_list=[-0.02, 0.32], imgname="j1302_awproject_iter2_tt0")
        self.mom8_creator(image=img1+'.image.tt1', range_list=[-0.07, 0.16], imgname="j1302_awproject_iter2_tt1")

        #################################################
        # %% Generate Images [test_j1302_awproject] end @
        #################################################

        # save results for future analysis
        np.save(self.id()+'.tt0stats.npy', tt0stats)
        np.save(self.id()+'.tt1stats.npy', tt1stats)
        np.save(self.id()+'.alphastats.npy', alphastats)
        np.save(self.id()+'.beamstats.npy', beamstats_curr)
        np.save(self.id()+'.tcleanrecs.npy', records)

        # (f) Runtimes not significantly different relative to previous runs
        # if ParallelTaskHelper.isMPIEnabled():
        #     # runtime with MPI -n 8
        #     successr, reportr = self.check_runtime(starttime, 9594)
        # else:
        #     successr, reportr = self.check_runtime(starttime, stats621['runtime'])

        # TODO start obeying on-axis once on-axis passes
        report  = "".join([report1, report3, report5, report6, report7])
        success = success1 and success3 and success5 and success6 and success7 and self.th.check_final(report)
        # report  = "".join([report0, report1, report2, report3, report4, report5, report6, report7])
        # success = success0 and success1 and success2 and success3 and success4 and success5 and success6 and success7 and self.th.check_final(report)

        add_to_dict(self, output = test_dict, dataset = "J1302-12fields.ms")
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = self.mom8_images
        test_dict[test_name]['taskcall'] = self.clean_taskcall(test_dict[test_name]['taskcall'], locals())

        self.assertTrue(success, msg=self.th.extract_failing_lines(report))

    # Test 3
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Tclean crashes with mpicasa+mosaic gridder+stokes imaging.")
    @stats_dict(test_dict)
    def test_j1302_mosaic_cube(self):
        """ [j1302] test_j1302_mosaic_cube """
        ######################################################################################
        # Should match values for "Cube" in the "Values to be compared"
        ######################################################################################

        ####################################################
        # %% Set local vars [test_j1302_mosaic_cube] start @
        ####################################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        test_name = self._testMethodName
        data_path_dir  = 'J1302/Stakeholder-test-mosaic-cube-data'
        basetname = 'J1302_mosaic_cube'
        masks = ['combined.mask', 'QLcatmask.mask']
        quick_masks = ['combined_1000.mask', 'QLcatmask_1000.mask']

        if not quick_test:
            imsize=4000
            self.prepData(self.vis, data_path_dir, *masks, partial_results_dirname="partial_results_test_j1302_mosaic_cube")
            for i in range(len(masks)):
                os.system(f"mv {masks[i]} {basetname}_{masks[i]}")
        else:
            imsize=1000
            self.prepData(self.vis, data_path_dir, *quick_masks, partial_results_dirname="partial_results_test_j1302_mosaic_cube")
            for i in range(len(masks)):
                os.system(f"mv {quick_masks[i]} {basetname}_{masks[i]}")
        self.teardown_files += [f"{basetname}_{x}" for x in masks]

        spw_chans = ''
        rms = {'0': 0.0005612289083201638, '1': 0.0004997134396517132, '2': 0.0008543526968930547} # per-spw noise floor as measured from a full-scale image run, Range:[100,100],[3900,1900]
        starttime = datetime.now()

        # reference frequence to use per spectral window (spw)
        # The datasets were reduced to take up less space. This included dropping the unused spws.
        # The original spws were '2', '8', and '14'. These got remapped to '0', '1', and '2', respectively.
        refFreqDict  = {
            '0' :  '2.028GHz',
            '1' :  '2.796GHz',
            '2' :  '3.564GHz'
        }

        ###################################################
        # %% Set local vars [test_j1302_mosaic_cube] end  @
        # %% Run tclean [test_j1302_mosaic_cube] start    @
        ###################################################

        def iname(image_iter, spw, stokes):
            #return 'J1302_'+image_iter+'_'+spw.replace('~','-')+'_'+stokes
            return basetname+'_'+image_iter+'_'+spw.replace('~','-')+'_'+stokes

        tstobj = self
        common_args = {'vis':'J1302-12fields.ms', 'uvrange':'<12km', 'intent':'OBSERVE_TARGET#UNSPECIFIED', \
                       'imsize':imsize, 'cell':'0.6arcsec', 'phasecenter':'13:03:13.874 -10.51.16.73', \
                       'gridder':'mosaic', 'conjbeams':False, 'mosweight':False, 
                       'usepointing':True, \
                       'pblimit':0.1, \
                       'deconvolver':'mtmfs', 'nterms':1, 'smallscalebias':0.4, \
                       'weighting':'briggs', 'robust':1.0,\
                       'threshold':0.0, 'cyclefactor':3.0, 'restart':True, 'interactive':False, 'parallel':self.parallel}

        def run_tclean(imagename=None, spw='', stokes='I', niter=None, datacolumn='corrected',reffreq=None,  mask='', 
              usemask='user', pbmask=0.0, nsigma=2.0, cycleniter=500, scales=[0],
              calcres=True, calcpsf=True, psfcutoff=0.5, savemodel='none', common_args=common_args):
            params = locals()
            finalparams = self.filter_runtclean_parameters(params, common_args)
            # add images to self.imgs for clean-up
            if ('imagename' in finalparams):
                img = finalparams['imagename']
                if (img not in self.imgs):
                    self.imgs.append(img)
            return tclean(**finalparams)

        spws         = [ '0','1','2']
        stokesParams = ['IQUV']
        spwstats     = { '0': {'freq': 2.028, 'IQUV': np.array([0.0,0.0,0.0,0.0]), 'beam': np.array([0.0,0.0,0.0])},
                         '1': {'freq': 2.796, 'IQUV': np.array([0.0,0.0,0.0,0.0]), 'beam': np.array([0.0,0.0,0.0])},
                         '2': {'freq': 3.594, 'IQUV': np.array([0.0,0.0,0.0,0.0]), 'beam': np.array([0.0,0.0,0.0])} }

        # These values are lifted from logs with John Tobin's original test scripts on the CAS-12427 ticket.
        # They are only used as a comparison, to make sure we use all the same exact parameter vlaues.
      #  script_pars_vals_0 = {
      #      '0': { 'IQUV': self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False) },
       #     '1': { 'IQUV': self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False) },
       #     '2': { 'IQUV': self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False) },
       # }
        #script_pars_vals_1 = {
        #    '0': { 'IQUV': self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
#
        #   '1': { 'IQUV': self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
        #    '2': { 'IQUV': self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
        #}
        #script_pars_vals_2 = {
        #    '0': { 'IQUV': self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
        #    '1': { 'IQUV': self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
        #    '2': { 'IQUV': self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
       #}
       # script_pars_vals_3 = {
       #     '0': { 'IQUV': self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
       #     '1': { 'IQUV': self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
       #     '2': { 'IQUV': self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
        #}

        records = {}
        for spw in spws:
            spw_str = f"{spw}{spw_chans}"
            if spw not in records:
                records[spw] = [0]*4
            r = records[spw]
            for stokes in stokesParams:
                # initialize iter2, no cleaning
                image_iter='iter2'
                imagename = iname(image_iter, spw, stokes)
                r[0] = run_tclean( imagename=imagename, datacolumn='corrected', niter=0, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw] )

                # resume iter2 with QL mask
                r[1] = run_tclean( imagename=imagename, datacolumn='corrected', scales=[0,5,12], nsigma=3.0, niter=20000, cycleniter=500,
                                   mask=f"{basetname}_QLcatmask.mask", calcres=False, calcpsf=False, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw] )

                # resume iter2 with combined mask
                os.system('rm -rf *.workdirectory')
                os.system('rm -rf *iter2*.mask')
                r[2] = run_tclean( imagename=imagename, datacolumn='corrected', scales=[0,5,12], nsigma=3.0, niter=20000, cycleniter=500,
                                   mask=f"{basetname}_combined.mask", calcres=False, calcpsf=False, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw])

                #os.system('rm -rf iter2*.mask')
                r[3] = run_tclean( imagename=imagename, datacolumn='corrected', scales=[0,5,12], nsigma=4.5, niter=20000, cycleniter=100,
                                   mask="", calcres=False, calcpsf=False, usemask='pb', pbmask=0.4, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw] )

                halfsize = round(imsize / 4000 * 2000)
                box      = f"{halfsize},{halfsize},{halfsize},{halfsize}"
                tt0statsI=imstat(imagename=imagename+'.image.tt0',box=box,stokes='I')
                tt0statsQ=imstat(imagename=imagename+'.image.tt0',box=box,stokes='Q')
                tt0statsU=imstat(imagename=imagename+'.image.tt0',box=box,stokes='U')
                tt0statsV=imstat(imagename=imagename+'.image.tt0',box=box,stokes='V')
                rmstt0statsI=imstat(imagename=imagename+'.residual.tt0',stokes='I')
                rmstt0statsQ=imstat(imagename=imagename+'.residual.tt0',stokes='Q')
                rmstt0statsU=imstat(imagename=imagename+'.residual.tt0',stokes='U')
                rmstt0statsV=imstat(imagename=imagename+'.residual.tt0',stokes='V')
                spwstats[spw]['IQUV']=np.squeeze(np.array([tt0statsI['max'],tt0statsQ['max'],tt0statsU['max'],tt0statsV['max']]))
                spwstats[spw]['rmsIQUV']=np.squeeze(np.array([rmstt0statsI['rms'],rmstt0statsQ['rms'],rmstt0statsU['rms'],rmstt0statsV['rms']]))
                spwstats[spw]['SNR']=((spwstats[spw]['IQUV']/spwstats[spw]['rmsIQUV'])**2)**0.5

                header = imhead(imagename+'.psf.tt0')
                beam   = header['perplanebeams']['beams']['*0']['*0']
                beamstats = np.squeeze(np.array([ beam['major']['value'], beam['minor']['value'], beam['positionangle']['value'] ]))
                spwstats[spw]['beam'] = beamstats

        ####################################################
        # %% Run tclean [test_j1302_mosaic_cube] end       @
        # %% Math stuff [test_j1302_mosaic_cube] start     @
        ####################################################

        spwlist=list(spwstats.keys())
        nspws=len(spwlist)
        freqs=np.zeros(nspws)
        fluxes=np.zeros(nspws)
        for i in range(nspws):
           freqs[i]=spwstats[spwlist[i]]['freq']
           fluxes[i]=spwstats[spwlist[i]]['IQUV'][0]

        logfreqs=np.log10(freqs)
        logfluxes=np.log10(fluxes)

        from scipy.optimize import curve_fit

        def func(x, a, b):
           nu_0=3.0
           return a*(x-np.log10(nu_0))+b

        popt, pcov = curve_fit(func, logfreqs, logfluxes)

        #############################################################
        # %% Math stuff [test_j1302_mosaic_cube] end                @
        # %% Compare Expected Values [test_j1302_mosaic_cube] start @
        #############################################################

        stats621 = {
            'full': {
                'F_nu': 0.31231569,
                'alpha': 0.03663584,
                #'runtime': 7167
            },
            '1/4': {
                'F_nu': 0.28358972,
                'alpha': -0.39547056,
                #'runtime': 902
            }
        }
        stats621 = stats621['full'] if (not quick_test) else stats621[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        # N/A: no pbcore, rms, or subim images are created for this test

        # (g) Fit F_nu0 and Alpha from three cube planes and compare: 6.2.1, on-axis
        # compare to alpha (ground truth), and 6.2.1 (fitted for spws 2, 8, 14 from mosaic gridder in CASA 6.2.1)
        alpha = popt[0]
        f_nu0 = 10**popt[1]
        curr_stats        = np.squeeze(np.array([f_nu0,            alpha])) # [flux density, alpha]
        onaxis_stats      = np.array([           0.3337,          -0.0476])
        casa621_stats     = np.array([           stats621['F_nu'], stats621['alpha']])
        success0, report0 = self.check_metrics_flux(curr_stats[0], onaxis_stats[0],  valname="Frac Diff F_nu vs. on-axis", rms_or_std=np.mean(list(rms.values())))
        success1, report1 = self.check_metrics_flux(curr_stats[0], casa621_stats[0], valname="Frac Diff F_nu vs. 6.2.1 image", rms_or_std=np.mean(list(rms.values())))
        success2, report2 = self.check_metrics_alpha_fitted(curr_stats[1], onaxis_stats[1], valname="Abs Diff alpha vs. on-axis", pcov=pcov)
        success3, report3 = self.check_metrics_alpha_fitted(curr_stats[1], casa621_stats[1], valname="Abs Diff alpha vs. 6.2.1 image", pcov=pcov)

        spwstats_621={
            'full': {
                '0': { 'freq': 2.028,
                       'IQUV':    np.array([ 0.30234554,     -0.00169584,     -0.00040815, -0.00173693]),
                       'beam':    np.array([ 4.62908077,     3.75791144,      8.6552124]),
                       'rmsIQUV': np.array([ 0.00058094,     0.000532,        0.00052686,  0.00055693]),
                       'SNR':     np.array([ 520.43911237,   3.18764363,      0.77467731,  3.11873003]) },
                '1': { 'freq': 2.796,
                       'IQUV':    np.array([ 0.324628860,    -1.01364225e-04, -2.73056852e-04, -9.81398509e-04]),
                       'beam':    np.array([ 3.31802654,     2.81110954,      11.9341383]),
                       'rmsIQUV': np.array([ 0.00060886,     0.00056572,      0.00056446,  0.00058757]),
                       'SNR':     np.array([ 5.33173526e+02, 0.179176112,     0.483746380, 1.67025396e+00]) },
                '2': { 'freq': 3.594,
                       'IQUV':    np.array([ 0.30719528,     0.00112637,      -0.00234151, -0.00072421]),
                       'beam':    np.array([ 1.99766636,     1.61410868,      8.42595863]),
                       'rmsIQUV': np.array([ 0.00120132,     0.00096404,      0.00095636,  0.00094095]),
                       'SNR':     np.array([ 255.71516847,   1.16838666,      2.44836242,  0.76965535]) }
            },
            '1/4': {
                '0': { 'freq': 2.028,
                       'IQUV':    np.array([ 0.32124895,     -0.00071501,    -0.00053517,    -0.00094193 ]),
                       'beam':    np.array([ 4.70300007,     3.75634837,     9.51373291 ]),
                       'rmsIQUV': np.array([ 0.00066714,     0.00051439,     0.00050774,     0.00054061 ]),
                       'SNR':     np.array([ 481.53417799,   1.39002728,     1.05402692,     1.74234655 ]) },
                '1': { 'freq': 2.796,
                       'IQUV':    np.array([ 0.312349528,    6.64699066e-04, 6.29942224e-05, -7.18787895e-04 ]),
                       'beam':    np.array([ 3.35121107,     2.79194474,     12.17673683 ]),
                       'rmsIQUV': np.array([ 0.00060982,     0.0004657,      0.00046498,     0.00047679 ]),
                       'SNR':     np.array([ 5.12196761e+02, 1.42730885,     0.135478162,    1.50756306 ]) },
                '2': { 'freq': 3.594,
                       'IQUV':    np.array([ 0.25404328,     0.00060382,     -0.00169632,    -0.00073403 ]),
                       'beam':    np.array([ 2.01960516,     1.63659084,     8.21656036 ]),
                       'rmsIQUV': np.array([ 0.000969,       0.00067771,     0.00066761,     0.00065421 ]),
                       'SNR':     np.array([ 262.17179635,   0.89097026,     2.54089228,     1.12200519 ])
                }
            }
        }
        spwstats_621 = spwstats_621['full'] if (not quick_test) else spwstats_621[f"1/4"]

        spwstats_onaxis={
          '0': { 'IQUV':    np.array([3.09755385e-01, -1.39351614e-04, -1.01510414e-04,  3.48565959e-06]),
                 'rmsIQUV': np.array([9.06101077e-05,  6.18512027e-05,  6.08855185e-05,  6.16381383e-05]),
                 'SNR':     np.array([3.41855222e+03,  2.25301381e+00,  1.66723411e+00,  5.65503710e-02]),
                 'beam':    np.array([4.39902544, 2.97761726, -1.9803896]),
                 'freq': 2.028},
          '1': { 'IQUV':    np.array([3.33031625e-01, -1.70329688e-04, -5.54503786e-05,  1.51384829e-05]),
                 'rmsIQUV': np.array([6.83661207e-05,  5.41465193e-05,  5.37294874e-05,  5.48078784e-05]),
                 'SNR':     np.array([4.87129621e+03,  3.14571813e+00,  1.03202880e+00,  2.76209979e-01]),
                 'beam':    np.array([3.17962098, 2.23836327, -3.84393311]),
                 'freq': 2.796},
          '2': { 'IQUV':    np.array([3.26020330e-01, -1.19368524e-04,  1.94136555e-05, -3.02872763e-06]),
                 'rmsIQUV': np.array([6.07943859e-05,  4.44997526e-05,  4.54793500e-05,  4.45020873e-05]),
                 'SNR':     np.array([5.36267166e+03,  2.68245365e+00,  4.26867480e-01,  6.80581028e-02]),
                 'beam':    np.array([2.4293716, 1.61365998, -3.72186279]),
                 'freq': 3.564}
        }

        success4 = []
        report4 = []
        for spw in spws:
            # (h) IQUV flux densities of all three spws: 6.2.1
            successN, reportN = self.check_metrics_flux(spwstats[spw]['IQUV'], spwstats_621[spw]['IQUV'],    valname=f"Stokes Comparison (spw {spw}), Frac Diff IQUV vs 6.2.1", rms_or_std=np.mean(list(rms.values())))
            success4.append(successN)
            report4.append(reportN)
            # (i) IQUV flux densities of all three spws: on-axis measurements
            successN, reportN = self.check_metrics_flux(spwstats[spw]['IQUV'], spwstats_onaxis[spw]['IQUV'], valname=f"Stokes Comparison (spw {spw}), Frac Diff IQUV vs on-axis", rms_or_std=np.mean(list(rms.values())))
            # TODO start obeying on-axis once on-axis passes
            # success4.append(successN)
            # report4.append(reportN)
            # (j) Beam of all three spws:                6.2.1
            successN, reportN = self.check_fracdiff(spwstats[spw]['beam'], spwstats_621[spw]['beam'],        valname=f"Stokes Comparison (spw {spw}), Frac Diff Maj, Min, PA vs 6.2.1")
            success4.append(successN)
            report4.append(reportN)

        ###########################################################
        # %% Compare Expected Values [test_j1302_mosaic_cube] end @
        # %% Generate Images [test_j1302_mosaic_cube] start       @
        ###########################################################

        #self.mom8_creator(image=iname('iter2', '0', 'IQUV')+'.image.tt0', range_list=[0.002, 0.32], imgname="j1302_mosaic_cube_spw2")
        #self.mom8_creator(image=iname('iter2', '1', 'IQUV')+'.image.tt0', range_list=[0.001, 0.32], imgname="j1302_mosaic_cube_spw8")
        #self.mom8_creator(image=iname('iter2', '2', 'IQUV')+'.image.tt0', range_list=[0.001, 0.30], imgname="j1302_mosaic_cube_spw14")
        self.mom8_creator(image=iname('iter2', '0', 'IQUV')+'.image.tt0', range_list=[-0.02, 0.32], imgname="j1302_mosaic_cube_spw2")
        self.mom8_creator(image=iname('iter2', '1', 'IQUV')+'.image.tt0', range_list=[-0.01, 0.32], imgname="j1302_mosaic_cube_spw8")
        self.mom8_creator(image=iname('iter2', '2', 'IQUV')+'.image.tt0', range_list=[-0.01, 0.30], imgname="j1302_mosaic_cube_spw14")

        ###################################################
        # %% Generate Images [test_j1302_mosaic_cube] end @
        ###################################################

        # save results for future analysis
        np.save(self.id()+'.fluxdens_alpha.npy', curr_stats)
        np.save(self.id()+'.freqs.npy', freqs)
        np.save(self.id()+'.fluxes.npy', fluxes)
        np.save(self.id()+'.spwstats.npy', spwstats)
        np.save(self.id()+'.tcleanrecs.npy', records)

        # (f) Runtimes not significantly different relative to previous runs
        # successr, reportr = self.check_runtime(starttime, stats621['runtime'])

        # TODO start obeying on-axis once on-axis passes
        report  = "".join([report1, report2, report3, *report4])
        success = success1 and success2 and success3 and all(success4) and self.th.check_final(report)
        # report  = "".join([report0, report1, report2, report3, *report4])
        # success = success0 and success1 and success2 and success3 and all(success4) and self.th.check_final(report)

        add_to_dict(self, output = test_dict, dataset = "J1302-12fields.ms")
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = self.mom8_images
        test_dict[test_name]['taskcall'] = self.clean_taskcall(test_dict[test_name]['taskcall'], locals())

        self.assertTrue(success, msg=self.th.extract_failing_lines(report))

    # Test 4
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Only run in serial, since John Tobin only executed this test in serial (see 01/12/22 comment on CAS-12427).")
    @stats_dict(test_dict)
    def test_j1302_ql(self):
        """ [j1302] test_j1302_ql """
        ######################################################################################
        # Should match values for "QL" in the "Values to be compared"
        ######################################################################################

        ###########################################
        # %% Set local vars [test_j1302_ql] start @
        ###########################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        test_name = self._testMethodName
        data_path_dir  = 'J1302/Stakeholder-test-mosaic-data'
        #img0 = 'VLASS1.2.ql.T08t20.J1302.10.2048.v1.I.iter0'
        #img1 = 'VLASS1.2.ql.T08t20.J1302.10.2048.v1.I.iter1'
        basetname = 'J1302_ql'
        img0 = f'{basetname}_iter0'
        img1 = f'{basetname}_iter1'
        self.prepData(self.vis, data_path_dir, partial_results_dirname="partial_results_test_j1302_ql")
        imsize = 7290
        rms = 0.00034846254286391285 # noise floor as measured from a full-scale image run, Range: [3000,3000],[6990,3600]

        starttime = datetime.now()
        if quick_test:
            imsize = 1822

        #########################################
        # %% Set local vars [test_j1302_ql] end @
        # %% Run tclean [test_j1302_ql] start   @
        #########################################

        records = []
        tstobj = self
        common_args = {'vis':'J1302-12fields.ms', 'uvrange':'<12km', 'intent':'OBSERVE_TARGET#UNSPECIFIED', \
                       'imsize':imsize, 'cell':'1.0arcsec', 'phasecenter':'13:03:13.874 -10.51.16.73', \
                       'datacolumn':'data',\
                       'spw':'', 'reffreq':'3.0GHz', \
                       'gridder':'mosaic', 'conjbeams':False, 'mosweight':False, 
                       'pblimit':0.02, \
                       'deconvolver':'mtmfs', 'nterms':2, 'smallscalebias':0.4, \
                       'weighting':'briggs', 'robust':1.0,\
                       'threshold':0.0, 'restart':True, 'interactive':False, 'parallel':self.parallel}

        def run_tclean(imagename=None, niter=None,  mask='', usemask='user', pbmask=0.0,  nsigma=0.0, cycleniter=-1, cyclefactor=1.0, scales=[0],
              calcres=True, calcpsf=True, savemodel='none', restoration=None, restoringbeam=None, common_args=common_args):
            params = locals()
            finalparams = self.filter_runtclean_parameters(params, common_args)
            # add images to self.imgs for clean-up
            if ('imagename' in finalparams):
                img = finalparams['imagename']
                if (img not in self.imgs):
                    self.imgs.append(img)
            return tclean(**finalparams)

        # These values are lifted from logs with John Tobin's original test scripts on the CAS-12427 ticket.
        # They are only used as a comparison, to make sure we use all the same exact parameter vlaues.
#        script_pars_vals_0 = self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='VLASS1.2.ql.T08t20.J1302.10.2048.v1.I.iter0', imsize=[7290, 7290], cell='1.0arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=False, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=360.0, pointingoffsetsigdev=[], pblimit=0.2, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.0, restoration=False, restoringbeam='common', pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[], niter=0, gain=0.1, threshold='0.0mJy', nsigma=0.0, cycleniter=-1, cyclefactor=1.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False)
#       script_pars_vals_1 = self.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='VLASS1.2.ql.T08t20.J1302.10.2048.v1.I.iter1', imsize=[7290, 7290], cell='1.0arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=False, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=360.0, pointingoffsetsigdev=[], pblimit=0.2, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.0, restoration=True, restoringbeam='common', pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=500, cyclefactor=2.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)

        run_tclean(imagename=img0, niter=0, restoration=False)
        if not use_partial_results:
            for ext in ['.weight.tt2', '.weight.tt0', '.psf.tt0', '.residual.tt0', '.weight.tt1', '.sumwt.tt2', '.psf.tt1', '.residual.tt1', '.psf.tt2', '.sumwt.tt1', '.model.tt0', '.pb.tt0', '.model.tt1', '.sumwt.tt0']:
                shutil.copytree(src=img0+ext, dst=img1+ext)
        run_tclean(imagename=img1, niter=20000, restoration=True, restoringbeam='common', nsigma=4.5, cycleniter=500, cyclefactor=2.0, calcres=False, calcpsf=False)

        if use_partial_results:
            os.system("rm -rf *.pbcor.tt0 *.subim")

        ###########################################
        # %% Run tclean [test_j1302_ql] end       @
        # %% Prepare Images [test_j1302_ql] start @
        ###########################################

        # hifv_pbcor(pipelinemode="automatic")
        for fromext,toext in [('.image.tt0','.image.pbcor.tt0'), ('.residual.tt0','.image.residual.pbcor.tt0')]:
            impbcor(imagename=img1+fromext, pbimage=img1+'.pb.tt0', outfile=img1+toext, mode='divide', cutoff=-1.0, stretch=False)
            self.check_img_exists(img1+toext)

        # hif_makermsimages(pipelinemode="automatic")
        imdev(imagename=img1+'.image.pbcor.tt0',
              outfile=img1+'.image.pbcor.tt0.rms',
              overwrite=True, stretch=False, grid=[10, 10], anchor='ref',
              xlength='60arcsec', ylength='60arcsec', interp='cubic', stattype='xmadm',
              statalg='chauvenet', zscore=-1, maxiter=-1)
        self.check_img_exists(img1+'.image.pbcor.tt0.rms')

        # hif_makecutoutimages(pipelinemode="automatic")
        blc = round(imsize / 7290 * 1785)
        urc = round(imsize / 7290 * 5506)
        for ext in ['.image.tt0', '.residual.tt0', '.image.pbcor.tt0', '.image.pbcor.tt0.rms', '.psf.tt0', '.image.residual.pbcor.tt0', '.pb.tt0']:
            imhead(imagename=img1+ext)
            imsubimage(imagename=img1+ext, outfile=img1+ext+'.subim', box=f"{blc},{blc},{urc},{urc}")
            self.check_img_exists(img1+ext+'.subim')

        ####################################################
        # %% Prepare Images [test_j1302_ql] end            @
        # %% Compare Expected Values [test_j1302_ql] start @
        ####################################################

        stats621 = {
            'full': {
                'F_nu': 0.3205489218235016,
                'beam': { 'maj': 3.12365127, 'min': 2.62010241, 'pos': 13.79291248 },
                #'runtime': 1805
            },
            '1/4': {
                'F_nu': 0.3184760808944702,
                'beam': { 'maj': 3.13571501, 'min': 2.6162672,  'pos': 13.37007141 },
                #'runtime': 383
            }
        }
        stats621 = stats621['full'] if (not quick_test) else stats621[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        success0, report0 = self.get_imgs_exist_results()

        # (a) tt0 vs 6.2.1, on-axis
        halfsize          = round(imsize / 7290 * 1860)
        imstat_vals       = imstat(imagename=img1+'.image.pbcor.tt0.subim',box=f"{halfsize},{halfsize},{halfsize},{halfsize}")
        curr_stats        = np.squeeze(np.array([imstat_vals['max']]))
        onaxis_stats      = np.array([           0.3337])
        casa621_stats     = np.array([           stats621['F_nu']])
        success1, report1 = self.check_metrics_flux(curr_stats, onaxis_stats, valname="Frac Diff F_nu vs. on-axis", rms_or_std=rms)
        success2, report2 = self.check_metrics_flux(curr_stats, casa621_stats, valname="Frac Diff F_nu vs. 6.2.1 image", rms_or_std=rms)

        # (b) tt1 vs 6.2.1, on-axis
        # no tt1 images for this test, skip

        # (c) alpha images
        # TODO
        # success3, report3 = ...

        # (d) beamsize comparison vs 6.2.1
        restbeam          = imhead(img1+'.image.pbcor.tt0.subim')['restoringbeam']
        beamstats_curr    = np.array([restbeam['major']['value'], restbeam['minor']['value'], restbeam['positionangle']['value']])
        beamstats_621     = np.array([stats621['beam']['maj'],    stats621['beam']['min'],    stats621['beam']['pos']])
        success4, report4 = self.check_fracdiff(beamstats_curr, beamstats_621, valname="Frac Diff Maj, Min, PA vs 6.2.1")

        ##################################################
        # %% Compare Expected Values [test_j1302_ql] end @
        # %% Generate Images [test_j1302_ql] start       @
        ##################################################

        #self.mom8_creator(image=img1+'.image.pbcor.tt0.subim', range_list=[0.001, 0.32], imgname="j1302_ql")
        self.mom8_creator(image=img1+'.image.pbcor.tt0.subim', range_list=[-0.1, 0.32], imgname="j1302_ql")

        ##########################################
        # %% Generate Images [test_j1302_ql] end @
        ##########################################

        # save results for future analysis
        np.save(self.id()+'.sourceflux.npy', curr_stats)
        np.save(self.id()+'.beamstats.npy', beamstats_curr)
        np.save(self.id()+'.tcleanrecs.npy', records)

        # (f) Runtimes not significantly different relative to previous runs
        # successr, reportr = self.check_runtime(starttime, stats621['runtime'])
        report  = "".join([report0, report1, report2, report4])
        success = success1 and success2 and success4 and self.th.check_final(report)

        add_to_dict(self, output = test_dict, dataset = "J1302-12fields.ms")
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = self.mom8_images
        test_dict[test_name]['taskcall'] = self.clean_taskcall(test_dict[test_name]['taskcall'], locals())

        self.assertTrue(success, msg=self.th.extract_failing_lines(report))


##############################################
##############################################
class test_j1927(StkUnitTest):

    def setUp(self):
        super().setUp()
        self.vis = 'J1927_12fields.ms'
        self.phasecenter = '19:27:30.443 +61.17.32.898'
        self.im = imager()
        self.parallel = False
        if ParallelTaskHelper.isMPIEnabled():
            self.parallel = True


    def tearDown(self):
        generate_weblog("tclean_VLASS_1v2_pipeline", test_dict)
        super().tearDown()
        if not cache_partial_results:
            self.delData() 

    # Test 5
    @stats_dict(test_dict)
    def test_j1927_mtmfs(self):
        """ [j1927] test_j1927_mtmfs """
        ######################################################################################
        # Should match values for "Stokes I" in the "Values to be compared"
        ######################################################################################
        # "noncube" to allow give this test a unique prefix, for running with runtest
        # For example: runtest.py -v test_vlass_1v2.py[test_j1927_mtmfs]

        ##############################################
        # %% Set local vars [test_j1927_mtmfs] start @
        ##############################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        test_name = self._testMethodName
        data_path_dir  = 'J1927/J1927-stakeholdertest-mosaic-data'
        #img0 = 'J1927_iter2'
        basetname = 'J1927_mtmfs'
        img0 = f'{basetname}_iter2'
        masks = ['secondmask.mask', 'QLcatmask.mask']
        quick_masks = ['secondmask_1000.mask', 'QLcatmask_1000.mask']

        if not quick_test:
            imsize=4000
            self.prepData(self.vis, data_path_dir, *masks, partial_results_dirname="partial_results_test_j1927_mtmfs")
            for i in range(len(masks)):
                os.system(f"mv {masks[i]} {basetname}_{masks[i]}")
        else:
            imsize=1000
            self.prepData(self.vis, data_path_dir, *quick_masks, partial_results_dirname="partial_results_test_j1927_mtmfs")
            for i in range(len(masks)):
                os.system(f"mv {quick_masks[i]} {basetname}_{masks[i]}")
        self.teardown_files += [f"{basetname}_{x}" for x in masks]

        spw = ''
        rms = [0.0001483304420688553, 0.0007968044018725578] # tt0, tt1 noise floor as measured from a full-scale image run, Range: [500,500],[3400,1900]
        starttime = datetime.now()

        #############################################
        # %% Set local vars [test_j1927_mtmfs] end  @
        # %% Prepare masks [test_j1927_mtmfs] start @
        #############################################

        # combine first and 2nd order masks
        if not use_partial_results:
            immath(imagename=[f'{basetname}_secondmask.mask',f'{basetname}_QLcatmask.mask'],expr='IM0+IM1',outfile=f'{basetname}_sum_of_masks.mask')
            self.im.mask(image=f'{basetname}_sum_of_masks.mask',mask=f'{basetname}_combined.mask',threshold=0.5)
            self.teardown_files += [f'{basetname}_sum_of_masks.mask', f'{basetname}_combined.mask']

        ###########################################
        # %% Prepare masks [test_j1927_mtmfs] end @
        # %% Run tclean [test_j1927_mtmfs] start  @
        ###########################################

        records = []
        tstobj = self
        common_args = {'vis':'J1927_12fields.ms', 'uvrange':'<12km', 'intent':'OBSERVE_TARGET#UNSPECIFIED', \
                       'imsize':imsize, 'cell':'0.6arcsec', 'phasecenter':'19:27:30.443 +61.17.32.898', \
                       'spw':'', 'reffreq':'3.0GHz', \
                       'gridder':'mosaic', 'mosweight':False,\
                       'pblimit':0.1, \
                       'deconvolver':'mtmfs', 'nterms':2, 'smallscalebias':0.4, \
                       'weighting':'briggs', 'robust':1.0,\
                       'threshold':0.0, 'cyclefactor':3.0, 'restart':True, 'interactive':False, 'parallel':self.parallel}

        def run_tclean(imagename=None, niter=None, datacolumn='corrected', mask='', usemask='user', pbmask=0.0, nsigma=2.0, cycleniter=500, scales=[0], calcres=False, calcpsf=False, savemodel='none', common_args=common_args):
            params = locals()
            finalparams = self.filter_runtclean_parameters(params, common_args)
            # add images to self.imgs for clean-up
            if ('imagename' in finalparams):
                img = finalparams['imagename']
                if (img not in self.imgs):
                    self.imgs.append(img)
            return tclean(**finalparams)

        # These values are lifted from logs with John Tobin's original test scripts on the CAS-12427 ticket.
        # They are only used as a comparison, to make sure we use all the same exact parameter vlaues.
#        script_pars_vals_0 = self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False)
#        script_pars_vals_1 = self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)
#       script_pars_vals_2 = self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='J1927_iter2', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='modelcolumn', calcres=False, calcpsf=False, parallel=False)
#       script_pars_vals_3 = self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)
#        script_pars_vals_4 = self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)

        # initialize iter2, no cleaning
        run_tclean(imagename=img0, niter=0, datacolumn='corrected', calcres=True, calcpsf=True )
        # # resume iter2 with QL mask
        run_tclean(imagename=img0, niter=20000, datacolumn='corrected', mask=f"{basetname}_QLcatmask.mask", nsigma=3.0, scales=[0,5,12] )

        # save model column, doesn't happen here in acutal VLASS pipeline, but makes sure functionality works.
        run_tclean(imagename=img0,  niter=0,     datacolumn='data',      savemodel='modelcolumn' )

        # resume iter2 with combined mask, remove old mask first, pass new mask as parameter
        os.system(f"rm -rf {img0}.mask")
        run_tclean(imagename=img0,  niter=20000, datacolumn='corrected', mask=f"{basetname}_combined.mask",  nsigma=3.0, scales=[0,5,12] )

        # resume iter2 with pbmask, removed old mask first then specify pbmask in resumption of tclean
        os.system(f"rm -rf {img0}.mask")
        run_tclean(imagename=img0, niter=20000, datacolumn='corrected', usemask='pb', mask="", pbmask=0.4,   nsigma=4.5, scales=[0,5,12], cycleniter=100 )

        #######################################################
        # %% Run tclean [test_j1927_mtmfs] end                @
        # %% Compare Expected Values [test_j1927_mtmfs] start @
        #######################################################

        stats621 = {
            'full': {
                'tt0': 0.88789159,
                'tt1': 0.41305166,
                'alpha': 0.46520507,
                'beam': { 'maj': 2.44883513, 'min': 2.04789495, 'pos': -23.64226532 },
                #'runtime': 4415
            },
            '1/4': {
                'tt0': 0.87826717,
                'tt1': 0.32428792,
                'alpha': 0.36923608,
                'beam': { 'maj': 2.35327077, 'min': 1.97257483, 'pos': -23.63233376 },
                #'runtime': 436
            }
        }
        stats621 = stats621['full'] if (not quick_test) else stats621[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        # N/A for this test

        halfsize          = round(imsize / 4000 * 2000)
        box               = f"{halfsize},{halfsize},{halfsize},{halfsize}"
        tt0stats   = imstat(imagename=img0+'.image.tt0', box=box)
        tt1stats   = imstat(imagename=img0+'.image.tt1', box=box)
        alphastats = imstat(imagename=img0+'.alpha',     box=box)
        curr_stats    = np.squeeze(np.array([ tt0stats['max'], tt1stats['max'], alphastats['max']]))
        onaxis_stats  = np.array([            0.9509,          0.3601,          0.3796])
        casa621_stats = np.array([            stats621['tt0'], stats621['tt1'], stats621['alpha']])

        # (a) tt0 vs 6.2.1, on-axis
        success0, report0 = self.check_metrics_flux(curr_stats[0], onaxis_stats[0],  valname="Frac Diff tt0 vs. on-axis", rms_or_std=rms[0])
        success1, report1 = self.check_metrics_flux(curr_stats[0], casa621_stats[0], valname="Frac Diff tt0 vs. 6.2.1 image", rms_or_std=rms[0])

        # (b) tt1 vs 6.2.1, on-axis
        success2, report2 = self.check_metrics_flux(curr_stats[1], onaxis_stats[1],  valname="Frac Diff tt1 vs. on-axis", rms_or_std=rms[1])
        success3, report3 = self.check_metrics_flux(curr_stats[1], casa621_stats[1], valname="Frac Diff tt1 vs. 6.2.1 image", rms_or_std=rms[1])

        # (c) alpha images
        success4, report4 = self.check_metrics_alpha(curr_stats[2], onaxis_stats[2],  onaxis_stats[0], onaxis_stats[1], valname="Abs Diff alpha vs. on-axis", rmss_or_stds=rms)
        success5, report5 = self.check_metrics_alpha(curr_stats[2], casa621_stats[2], casa621_stats[0], casa621_stats[1], valname="Abs Diff alpha vs. 6.2.1 image", rmss_or_stds=rms)

        # (d) beamsize comparison vs 6.2.1
        restbeam          = imhead(img0+'.image.tt0')['restoringbeam']
        beamstats_curr    = np.array([restbeam['major']['value'], restbeam['minor']['value'], restbeam['positionangle']['value']])
        beamstats_621     = np.array([stats621['beam']['maj'],    stats621['beam']['min'],    stats621['beam']['pos']])
        success6, report6 = self.check_fracdiff(beamstats_curr, beamstats_621, valname="Frac Diff Maj, Min, PA vs 6.2.1")

        # (e) Confirm presence of model column in resultant MS
        success7, report7 = self.check_column_exists("MODEL_DATA")

        #####################################################
        # %% Compare Expected Values [test_j1927_mtmfs] end @
        # %% Generate Images [test_j1927_mtmfs] start       @
        #####################################################

        #self.mom8_creator(image=img0+'.image.tt0', range_list=[0, 0.88], imgname="j1927_mtmfs_tt0")
        #self.mom8_creator(image=img0+'.image.tt0', range_list=[0, 0.88], imgname="j1927_mtmfs_tt0")
        self.mom8_creator(image=img0+'.image.tt0', range_list=[-0.05, 0.88], imgname="j1927_mtmfs_tt0")
        self.mom8_creator(image=img0+'.image.tt1', range_list=[-0.05, 0.41], imgname="j1927_mtmfs_tt1")

        #############################################
        # %% Generate Images [test_j1927_mtmfs] end @
        #############################################

        # save results for future analysis
        np.save(self.id()+'.tt0tt1alpha.npy', curr_stats)
        np.save(self.id()+'.beamstats.npy', beamstats_curr)
        np.save(self.id()+'.tcleanrecs.npy', records)

        # (f) Runtimes not significantly different relative to previous runs
        # if ParallelTaskHelper.isMPIEnabled():
        #     # runtime with MPI -n 8
        #     successr, reportr = self.check_runtime(starttime, 7129)
        # else:
        #     successr, reportr = self.check_runtime(starttime, stats621['runtime'])

        # TODO start obeying on-axis once on-axis passes
        report  = "".join([report0, report1, report3, report4, report5, report6, report7])
        success = success0 and success1 and success3 and success4 and success5 and success6 and success7 and self.th.check_final(report)
        # report  = "".join([report0, report1, report2, report3, report4, report5, report6, report7, reportr])
        # success = success0 and success1 and success2 and success3 and success4 and success5 and success6 and success7 and successr and self.th.check_final(report)

        add_to_dict(self, output = test_dict, dataset = "J1927-12fields.ms")
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = self.mom8_images
        test_dict[test_name]['taskcall'] = self.clean_taskcall(test_dict[test_name]['taskcall'], locals())

        self.assertTrue(success, msg=self.th.extract_failing_lines(report))

    # N/A not implemented
    # def test_j1927_awproject(self):
    #     pass

    # Test 6
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Tclean crashes with mpicasa+mosaic gridder+stokes imaging.")
    @stats_dict(test_dict)
    def test_j1927_mosaic_cube(self):
        """ [j1927] test_j1927_mosaic_cube """
        ######################################################################################
        # Should match values for "Cube" in the "Values to be compared"
        ######################################################################################

        ####################################################
        # %% Set local vars [test_j1927_mosaic_cube] start @
        ####################################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        test_name = self._testMethodName
        data_path_dir  = 'J1927/J1927-stakeholdertest-mosaic-cube-data'
        basetname = 'J1927_mosaic_cube'
        masks = ['combined.mask', 'QLcatmask.mask']
        quick_masks = ['combined_1000.mask', 'QLcatmask_1000.mask']

        if not quick_test:
            imsize=4000
            self.prepData(self.vis, data_path_dir, *masks, partial_results_dirname="partial_results_test_j1927_mosaic_cube")
            for i in range(len(masks)):
                os.system(f"mv {masks[i]} {basetname}_{masks[i]}")
        else:
            imsize=1000
            self.prepData(self.vis, data_path_dir, *quick_masks, partial_results_dirname="partial_results_test_j1927_mosaic_cube")
            for i in range(len(masks)):
                os.system(f"mv {quick_masks[i]} {basetname}_{masks[i]}")
        self.teardown_files += [f'{basetname}_{x}' for x in masks]

        # rundir = "/users/bbean/dev/CAS-12427/src/casalith/build-casalith/work/linux/test_vlass_j1927_cube_unittest"
        # os.system(f"mv {rundir}/run_results/VLASS* {rundir}/nosedir/test_vlass_1v2/")
        spw_chans = ''
        rms = {'0': 0.0004637447465635465, '1': 0.0005062325702676312, '2': 0.0004403419438565781} # per-spw noise floor as measured from a full-scale image run, Range:[100,100],[3900,1900]
        starttime = datetime.now()

        # reference frequence to use per spectral window (spw)
        # The datasets were reduced to take up less space. This included dropping the unused spws.
        # The original spws were '2', '8', and '14'. These got remapped to '0', '1', and '2', respectively.
        refFreqDict  = {
            '0' :  '2.028GHz',
            '1' :  '2.796GHz',
            '2' :  '3.564GHz'
        }

        ###################################################
        # %% Set local vars [test_j1927_mosaic_cube] end  @
        # %% Prepare masks [test_j1927_mosaic_cube] start @
        ###################################################

        #################################################
        # %% Prepare masks [test_j1927_mosaic_cube] end @
        # %% Run tclean [test_j1927_mosaic_cube] start  @
        #################################################

        def iname(image_iter, spw, stokes):
        #    return 'J1927_'+image_iter+'_'+spw.replace('~','-')+'_'+stokes
            return basetname+'_'+image_iter+'_'+spw.replace('~','-')+'_'+stokes

        tstobj = self
        common_args = {'vis':'J1927_12fields.ms', 'uvrange':'<12km', 'intent':'OBSERVE_TARGET#UNSPECIFIED', \
                       'imsize':imsize, 'cell':'0.6arcsec', 'phasecenter':'19:27:30.443 +61.17.32.898', \
                       'gridder':'mosaic', 'mosweight':False, \
                       'pblimit':0.1, \
                       'deconvolver':'mtmfs', 'nterms':1, 'smallscalebias':0.4, \
                       'weighting':'briggs', 'robust':1.0,\
                       'threshold':0.0, 'cyclefactor':3.0, 'restart':True, 'interactive':False, 'parallel':self.parallel}

        def run_tclean(imagename=None, niter=None, spw='', stokes='I', datacolumn='corrected', reffreq=None, mask='', usemask='user', pbmask=0.0, cfcache='', wbawp=True,  nsigma=2.0, cycleniter=5000, scales=[0],
              calcres=True, calcpsf=True, psfcutoff=0.5, savemodel='none', common_args=common_args):
            params = locals()
            finalparams = self.filter_runtclean_parameters(params, common_args)
            # add images to self.imgs for clean-up
            if ('imagename' in finalparams):
                img = finalparams['imagename']
                if (img not in self.imgs):
                    self.imgs.append(img)
            return tclean(**finalparams)


        spws         = [ '0','1','2']
        stokesParams = ['IQUV']
        spwstats     = { '0': {'freq': 2.028, 'IQUV': np.array([0.0,0.0,0.0,0.0]), 'beam': np.array([0.0,0.0,0.0])},
                         '1': {'freq': 2.796, 'IQUV': np.array([0.0,0.0,0.0,0.0]), 'beam': np.array([0.0,0.0,0.0])},
                         '2': {'freq': 3.594, 'IQUV': np.array([0.0,0.0,0.0,0.0]), 'beam': np.array([0.0,0.0,0.0])} }

        # These values are lifted from logs with John Tobin's original test scripts on the CAS-12427 ticket.
        # They are only used as a comparison, to make sure we use all the same exact parameter vlaues.
        script_pars_vals_0 = {
            '0': { 'IQUV': self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False) },
            '1': { 'IQUV': self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False) },
            '2': { 'IQUV': self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False) },
        }
        script_pars_vals_1 = {
            '0': { 'IQUV': self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '1': { 'IQUV': self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '2': { 'IQUV': self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
        }
        script_pars_vals_2 = {
            '0': { 'IQUV': self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '1': { 'IQUV': self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '2': { 'IQUV': self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
        }
        script_pars_vals_3 = {
            '0': { 'IQUV': self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '1': { 'IQUV': self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '2': { 'IQUV': self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
        }

        records = {}
        for spw in spws:
            spw_str = f"{spw}{spw_chans}"
            if spw not in records:
                records[spw] = [0]*4
            r = records[spw]
            for stokes in stokesParams:
                # initialize iter2, no cleaning
                image_iter='iter2'
                imagename = iname(image_iter, spw, stokes)
                r[0] = run_tclean( imagename=imagename, datacolumn='corrected', niter=0, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw] )

                # # resume iter2 with QL mask
                r[1] = run_tclean( imagename=imagename, datacolumn='corrected', scales=[0,5,12], nsigma=3.0, niter=20000, cycleniter=500,   
                                   mask=f"{basetname}_QLcatmask.mask", calcres=False, calcpsf=False, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw] )

                # resume iter2 with combined mask
                os.system('rm -rf *.workdirectory')
                os.system('rm -rf *iter2*.mask')
                r[2] = run_tclean( imagename=imagename, datacolumn='corrected', scales=[0,5,12], nsigma=3.0, niter=20000, cycleniter=500,
                                   mask=f"{basetname}_combined.mask", calcres=False, calcpsf=False, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw])

                # os.system('rm -rf iter2*.mask')
                r[3] = run_tclean( imagename=imagename, datacolumn='corrected', scales=[0,5,12], nsigma=4.5, niter=20000, cycleniter=100,
                                   mask="", calcres=False, calcpsf=False, usemask='pb', pbmask=0.4, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw] )

                halfsize = round(imsize / 4000 * 2000)
                box      = f"{halfsize},{halfsize},{halfsize},{halfsize}"
                tt0statsI=imstat(imagename=imagename+'.image.tt0',box=box,stokes='I')
                tt0statsQ=imstat(imagename=imagename+'.image.tt0',box=box,stokes='Q')
                tt0statsU=imstat(imagename=imagename+'.image.tt0',box=box,stokes='U')
                tt0statsV=imstat(imagename=imagename+'.image.tt0',box=box,stokes='V')
                spwstats[spw]['IQUV']=np.squeeze(np.array([tt0statsI['max'],tt0statsQ['max'],tt0statsU['max'],tt0statsV['max']]))
                header=imhead(imagename+'.psf.tt0')['perplanebeams']['beams']['*0']['*0']
                beamstats=np.squeeze(np.array([header['major']['value'],header['minor']['value'],header['positionangle']['value']]))
                spwstats[spw]['beam']=beamstats

        ####################################################
        # %% Run tclean [test_j1927_mosaic_cube] end       @
        # %% Math stuff [test_j1927_mosaic_cube] start     @
        ####################################################

        spwlist=list(spwstats.keys())
        nspws=len(spwlist)
        freqs=np.zeros(nspws)
        fluxes=np.zeros(nspws)
        for i in range(nspws):
           freqs[i]=spwstats[spwlist[i]]['freq']
           fluxes[i]=spwstats[spwlist[i]]['IQUV'][0]

        logfreqs=np.log10(freqs)
        logfluxes=np.log10(fluxes)

        from scipy.optimize import curve_fit

        def func(x, a, b):
            nu_0=3.0
            return a*(x-np.log10(nu_0))+b

        popt, pcov = curve_fit(func, logfreqs, logfluxes)

        #############################################################
        # %% Math stuff [test_j1927_mosaic_cube] end                @
        # %% Compare Expected Values [test_j1927_mosaic_cube] start @
        #############################################################

        stats621 = {
            'full': {
                'F_nu': 0.88835657,
                'alpha': 0.41202253,
                #'runtime': 15329
            },
            '1/4': {
                'F_nu': 0.87965548,
                'alpha': 0.29201263,
                #'runtime': 1022
            }
        }
        stats621 = stats621['full'] if (not quick_test) else stats621[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        # N/A: no pbcore, rms, or subim images are created for this test

        # (g) Fit F_nu0 and Alpha from three cube planes and compare: 6.2.1, on-axis
        # compare to alpha (ground truth), and 6.2.1 (fitted for spws 2, 8, 14 from mosaic gridder in CASA 6.2.1)
        alpha = popt[0]
        f_nu0 = 10**popt[1]
        curr_stats        = np.squeeze(np.array([f_nu0,            alpha])) # [flux density, alpha]
        onaxis_stats      = np.array([           0.9509,           0.3601])
        casa621_stats     = np.array([           stats621['F_nu'], stats621['alpha']])
        success0, report0 = self.check_metrics_flux(curr_stats[0], onaxis_stats[0],  valname="Frac Diff F_nu vs. on-axis", rms_or_std=np.mean(list(rms.values())))
        success1, report1 = self.check_metrics_flux(curr_stats[0], casa621_stats[0], valname="Frac Diff F_nu vs. 6.2.1 image", rms_or_std=np.mean(list(rms.values())))
        success2, report2 = self.check_metrics_alpha_fitted(curr_stats[1], onaxis_stats[1],  valname="Abs Diff alpha vs. on-axis", pcov=pcov)
        success3, report3 = self.check_metrics_alpha_fitted(curr_stats[1], casa621_stats[1], valname="Abs Diff alpha vs. 6.2.1 image", pcov=pcov)

        spwstats_621= {
            'full': {
                '0': { 'freq': 2.028,
                       'IQUV': np.array([ 0.75540608, 0.00575879,  0.0007934 , -0.00628786]),
                       'beam': np.array([ 3.65377903, 3.03121424, -22.77662277]) },
                '1': { 'freq': 2.796,
                       'IQUV': np.array([ 0.86449653, 0.00522652, -0.00604673, -0.00736201]),
                       'beam': np.array([ 2.65529013, 2.2384088 , -26.70895958]) },
                '2': { 'freq': 3.594,
                       'IQUV': np.array([ 0.95604223, 0.00543462, -0.00965639, -0.00357582]),
                       'beam': np.array([ 2.06172824, 1.69355631, -25.49251938]) }
            },
            '1/4': {
                '0': { 'freq': 2.028,
                       'IQUV': np.array([ 0.779013932, 6.74977247e-03, 4.55768779e-04, 5.20252250e-03 ]),
                       'beam': np.array([ 3.65848589,  3.01737332,     -24.03954887 ]) },
                '1': { 'freq': 2.796,
                       'IQUV': np.array([ 0.87592876,  0.00488591,     -0.00589986,    0.00847113 ]),
                       'beam': np.array([ 2.63813567,  2.21187878,     -26.08016586 ]) },
                '2': { 'freq': 3.594,
                       'IQUV': np.array([ 0.91885197,  0.00677892,     -0.01262855,    0.01252271 ]),
                       'beam': np.array([ 2.05659461,  1.67871821,     -25.18776894 ])

                }
            }
        }

        spwstats_621 = spwstats_621['full'] if (not quick_test) else spwstats_621[f"1/4"]

        success4 = []
        report4 = []
        for spw in spws:
            # (h) IQUV flux densities of all three spws:              6.2.1
            successN, reportN = self.check_metrics_flux(spwstats[spw]['IQUV'], spwstats_621[spw]['IQUV'], valname=f"Stokes Comparison (spw {spw}), Frac Diff IQUV vs 6.2.1", rms_or_std=np.mean(list(rms.values())))
            success4.append(successN)
            report4.append(reportN)
            # (i) IQUV flux densities of all three spws:              on-axis measurements
            # N/A: no no-axis measurements available in VLASS_mosaic_cube_stakeholder_test_script.py
            # (j) Beam of all three spws:                             6.2.1
            successN, reportN = self.check_fracdiff(spwstats[spw]['beam'], spwstats_621[spw]['beam'],     valname=f"Stokes Comparison (spw {spw}), Frac Diff Maj, Min, PA vs 6.2.1")
            success4.append(successN)
            report4.append(reportN)

        ###########################################################
        # %% Compare Expected Values [test_j1927_mosaic_cube] end @
        # %% Generate Images [test_j1927_mosaic_cube] start       @
        ###########################################################

        #self.mom8_creator(image=iname('iter2', '0', 'IQUV')+'.image.tt0', range_list=[0, 0.78], imgname="j1927_mosaic_cube_spw2")
        #self.mom8_creator(image=iname('iter2', '1', 'IQUV')+'.image.tt0', range_list=[0, 0.88], imgname="j1927_mosaic_cube_spw8")
        #self.mom8_creator(image=iname('iter2', '2', 'IQUV')+'.image.tt0', range_list=[0, 0.92], imgname="j1927_mosaic_cube_spw14")
        self.mom8_creator(image=iname('iter2', '0', 'IQUV')+'.image.tt0', range_list=[-0.05, 0.78], imgname="j1927_mosaic_cube_spw2")
        self.mom8_creator(image=iname('iter2', '1', 'IQUV')+'.image.tt0', range_list=[-0.05, 0.88], imgname="j1927_mosaic_cube_spw8")
        self.mom8_creator(image=iname('iter2', '2', 'IQUV')+'.image.tt0', range_list=[-0.05, 0.92], imgname="j1927_mosaic_cube_spw14")

        ###################################################
        # %% Generate Images [test_j1927_mosaic_cube] end @
        ###################################################

        # save results for future analysis
        np.save(self.id()+'.fluxdens_alpha.npy', curr_stats)
        np.save(self.id()+'.freqs.npy', freqs)
        np.save(self.id()+'.fluxes.npy', fluxes)
        np.save(self.id()+'.spwstats.npy', spwstats)
        np.save(self.id()+'.tcleanrecs.npy', records)

        # (f) Runtimes not significantly different relative to previous runs
        # successr, reportr = self.check_runtime(starttime, stats621['runtime'])
        report  = "".join([report0, report1, report2, report3, *report4])
        success = success0 and success1 and success2 and success3 and all(success4) and self.th.check_final(report)

        add_to_dict(self, output = test_dict, dataset = "J1927-12fields.ms")
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = self.mom8_images
        test_dict[test_name]['taskcall'] = self.clean_taskcall(test_dict[test_name]['taskcall'], locals())

        self.assertTrue(success, msg=self.th.extract_failing_lines(report))

    # Test 7
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Only run in serial, since John Tobin only executed this test in serial (see 01/12/22 comment on CAS-12427).")
    @stats_dict(test_dict)
    def test_j1927_ql(self):
        """ [j1927] test_j1927_ql """
        ######################################################################################
        # Should match values for "QL" in the "Values to be compared"
        ######################################################################################

        ###########################################
        # %% Set local vars [test_j1927_ql] start @
        ###########################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        test_name = self._testMethodName
        data_path_dir  = 'J1927/J1927-stakeholdertest-mosaic-data'
        basetname = 'J1927_ql'
        #img0 = 'VLASS1.2.ql.T26t15.J1927.10.2048.v1.I.iter0'
        #img1 = 'VLASS1.2.ql.T26t15.J1927.10.2048.v1.I.iter1'
        img0 = f'{basetname}_iter0'
        img1 = f'{basetname}_iter1'
        self.prepData(self.vis, data_path_dir, partial_results_dirname="partial_results_test_j1927_ql")
        # rundir = "/users/bbean/dev/CAS-12427/src/casalith/build-casalith/work/linux/test_vlass_j1927_QL_unittest"
        # os.system(f"mv {rundir}/run_results/VLASS* {rundir}/nosedir/test_vlass_1v2/")
        imsize = 7290
        rms = 0.00023228885125825126 # noise floor as measured from a full-scale image run, Range: [2800,2900],[4300,3600]

        starttime = datetime.now()
        if quick_test:
            imsize = 1822

        #########################################
        # %% Set local vars [test_j1927_ql] end @
        # %% Run tclean [test_j1927_ql] start   @
        #########################################

        records = []
        tstobj = self
        common_args = {'vis':'J1927_12fields.ms', 'uvrange':'<12km', 'intent':'OBSERVE_TARGET#UNSPECIFIED', \
                       'imsize':imsize, 'cell':'1.0arcsec', 'phasecenter':'19:27:30.443 +61.17.32.898', \
                       'reffreq':'3.0GHz', \
                       'datacolumn':'data',\
                       'gridder':'mosaic', 'conjbeams':False, 'mosweight':False, \
                       'restoringbeam':'common',\
                       'pblimit':0.2, \
                       'deconvolver':'mtmfs', 'nterms':2, 'smallscalebias':0.0, \
                       'weighting':'briggs', 'robust':1.0,\
                       'threshold':0.0, 'restart':True, 'interactive':False, 'parallel':self.parallel}

        def run_tclean(imagename=None, niter=None,  mask='', usemask='user', pbmask=0.0, nsigma=0.0, cycleniter=-1, cyclefactor=1.0, scales=[0],
              calcres=True, calcpsf=True, restoration=None, common_args=common_args):
            params = locals()
            finalparams = self.filter_runtclean_parameters(params, common_args)
            # add images to self.imgs for clean-up
            if ('imagename' in finalparams):
                img = finalparams['imagename']
                if (img not in self.imgs):
                    self.imgs.append(img)
            return tclean(**finalparams)

        # These values are lifted from logs with John Tobin's original test scripts on the CAS-12427 ticket.
        # They are only used as a comparison, to make sure we use all the same exact parameter vlaues.
        script_pars_vals_0 = self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='VLASS1.2.ql.T26t15.J1927.10.2048.v1.I.iter0', imsize=[7290, 7290], cell='1.0arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=False, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=360.0, pointingoffsetsigdev=[], pblimit=0.2, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.0, restoration=False, restoringbeam='common', pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[], niter=0, gain=0.1, threshold='0.0mJy', nsigma=0.0, cycleniter=-1, cyclefactor=1.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False)
        script_pars_vals_1 = self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='VLASS1.2.ql.T26t15.J1927.10.2048.v1.I.iter1', imsize=[7290, 7290], cell='1.0arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=False, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=360.0, pointingoffsetsigdev=[], pblimit=0.2, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.0, restoration=True, restoringbeam='common', pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=500, cyclefactor=2.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)

        run_tclean(imagename=img0, niter=0, restoration=False)
        if not use_partial_results:
            for ext in ['.weight.tt2', '.weight.tt0', '.psf.tt0', '.residual.tt0', '.weight.tt1', '.sumwt.tt2', '.psf.tt1', '.residual.tt1', '.psf.tt2', '.sumwt.tt1', '.model.tt0', '.pb.tt0', '.model.tt1', '.sumwt.tt0']:
                shutil.copytree(src=img0+ext, dst=img1+ext)
        run_tclean(imagename=img1, niter=20000, restoration=True,  nsigma=4.5, cycleniter=500, cyclefactor=2.0, calcres=False, calcpsf=False)

        if use_partial_results:
            os.system("rm -rf *.pbcor.tt0* *.subim")

        ###########################################
        # %% Run tclean [test_j1927_ql] end       @
        # %% Prepare Images [test_j1927_ql] start @
        ###########################################

        # hifv_pbcor(pipelinemode="automatic")
        for fromext,toext in [('.image.tt0','.image.pbcor.tt0'), ('.residual.tt0','.image.residual.pbcor.tt0')]:
            impbcor(imagename=img1+fromext, pbimage=img1+'.pb.tt0', outfile=img1+toext, mode='divide', cutoff=-1.0, stretch=False)
            self.check_img_exists(img1+toext)

        # hif_makermsimages(pipelinemode="automatic")
        imdev(imagename=img1+'.image.pbcor.tt0',
              outfile=img1+'.image.pbcor.tt0.rms',
              overwrite=True, stretch=False, grid=[10, 10], anchor='ref',
              xlength='60arcsec', ylength='60arcsec', interp='cubic', stattype='xmadm',
              statalg='chauvenet', zscore=-1, maxiter=-1)
        self.check_img_exists(img1+'.image.pbcor.tt0.rms')

        # hif_makecutoutimages(pipelinemode="automatic")
        blc = round(imsize / 7290 * 1785)
        urc = round(imsize / 7290 * 5506)
        for ext in ['.image.tt0', '.residual.tt0', '.image.pbcor.tt0', '.image.pbcor.tt0.rms', '.psf.tt0', '.image.residual.pbcor.tt0', '.pb.tt0']:
            imhead(imagename=img1+ext)
            imsubimage(imagename=img1+ext, outfile=img1+ext+'.subim', box=f"{blc},{blc},{urc},{urc}")
            self.check_img_exists(img1+ext+'.subim')

        ####################################################
        # %% Prepare Images [test_j1927_ql] end            @
        # %% Compare Expected Values [test_j1927_ql] start @
        ####################################################

        stats621 = {
            'full': {
                'F_nu': [0.9051663279533386],
                'beam': { 'maj': 2.45210147, 'min': 2.07915711, 'pos': -18.88167 },
                #'runtime': 4277
            },
            '1/4': {
                'F_nu': 0.8975918889045715,
                'beam': { 'maj': 2.4284029,  'min': 2.06255174, 'pos': -19.85708046 },
                #'runtime': 826
            }
        }
        stats621 = stats621['full'] if (not quick_test) else stats621[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        success0, report0 = self.get_imgs_exist_results()

        # (a) tt0 vs 6.2.1, on-axis
        halfsize          = round(imsize / 7290 * 1860)
        box               = f"{halfsize},{halfsize},{halfsize},{halfsize}"
        imstat_vals       = imstat(imagename=img1+'.image.pbcor.tt0.subim', box=box)
        curr_stats        = np.squeeze(np.array([imstat_vals['max']]))
        onaxis_stats      = np.array([           0.9509])
        casa621_stats     = np.array([           stats621['F_nu']])
        success1, report1 = self.check_metrics_flux(curr_stats, onaxis_stats,  valname="Frac Diff F_nu vs. on-axis", rms_or_std=rms)
        success2, report2 = self.check_metrics_flux(curr_stats, casa621_stats, valname="Frac Diff F_nu vs. 6.2.1 image", rms_or_std=rms)

        # (b) tt1 vs 6.2.1, on-axis
        # no tt1 images for this test, skip

        # (c) alpha images
        # TODO
        # success3, report3 = ...

        # (d) beamsize comparison vs 6.2.1
        restbeam          = imhead(img1+'.image.pbcor.tt0.subim')['restoringbeam']
        beamstats_curr    = np.array([restbeam['major']['value'], restbeam['minor']['value'], restbeam['positionangle']['value']])
        beamstats_621     = np.array([stats621['beam']['maj'],    stats621['beam']['min'],    stats621['beam']['pos']])
        success4, report4 = self.check_fracdiff(beamstats_curr, beamstats_621, valname="Frac Diff Maj, Min, PA vs 6.2.1")

        ##################################################
        # %% Compare Expected Values [test_j1927_ql] end @
        # %% Generate Images [test_j1927_ql] start       @
        ##################################################

        #self.mom8_creator(image=img1+'.image.pbcor.tt0.subim', range_list=[0, 0.90], imgname="j1927_ql")
        self.mom8_creator(image=img1+'.image.pbcor.tt0.subim', range_list=[-0.1, 0.90], imgname="j1927_ql")

        ##########################################
        # %% Generate Images [test_j1927_ql] end @
        ##########################################

        # save results for future analysis
        np.save(self.id()+'.sourceflux.npy', curr_stats)
        np.save(self.id()+'.beamstats.npy', beamstats_curr)

        # (f) Runtimes not significantly different relative to previous runs
        #successr, reportr = self.check_runtime(starttime, stats621['runtime'])
        report  = "".join([report0, report1, report2, report4])
        success = success1 and success2 and success4 and self.th.check_final(report)

        add_to_dict(self, output = test_dict, dataset = "J1927-12fields.ms")
        test_dict[test_name]['report'] = report
        test_dict[test_name]['images'] = self.mom8_images
        test_dict[test_name]['taskcall'] = self.clean_taskcall(test_dict[test_name]['taskcall'], locals())

        self.assertTrue(success, msg=self.th.extract_failing_lines(report))

##############################################
##############################################

if __name__ == '__main__':
    unittest.main()
