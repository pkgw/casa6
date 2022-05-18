##########################################################################
#
# Run the tests as described in https://open-confluence.nrao.edu/pages/viewpage.action?spaceKey=CASA&title=Requirements+for+VLASS+Imaging+Pipeline+Stakeholders+Tests
# There wasn't a great definition of the work to be done, so much of it is being interpretted by me (BGB).
#
##########################################################################
#
# The basic idea:
#  Compare values from the most recent build of casa to other known values. Check for differences.
#
# Values to compare current versions against:
#  1. on-axis      ---   "ground truth", I think from the (fluxscale?) task
#  2. CASA 6.1.3   ---   pipeline approved version of casa
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
#   a. tt0:                                                    6.1.3, on-axis
#   b. tt1:                                                    6.1.3, on-axis
#   c. alpha:                                                  6.1.3, on-axis
#   d. beamsize comparison:                                    6.1.3
#  Stokes I:
#   e. Confirm presence of model column in resultant MS
#  Stokes I and Cube:
#   f. Runtimes not significantly different relative to previous runs
#  Cube:
#   g. Fit F_nu0 and Alpha from three cube planes and compare: 6.1.3, on-axis
#   h. IQUV flux densities of all three spws:                  6.1.3
#   i. IQUV flux densities of all three spws:                  on-axis measurements
#   j. Beam of all three spws:                                 6.1.3
#  QL:
#   k. flux density of Calibrator source:                      6.1.3, on-axis
#  All:
#   l. Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
# Note on (i): on-axis values derived from images of calibrator using standard gridder of pointed data on calibrator
#              ^-- does this mean these are also "ground truth" values? BGB211222
#
# Other goals:
#  AUX1: try to make these tests compatible with jupyter, although this requirement may be dropped due to the time it takes to run these tests
#   https://github.com/casangi/stakeholder/
#   https://colab.research.google.com/drive/1IdnKmhNWonX1r1zJ5oCTExVMSMgTHjOm?usp=sharing#scrollTo=J88MWrnHFXZP
#   https://docs.google.com/document/d/1rKktZdg4IwV5Nr_giVRcFlvkFkZ2NxTv-BO6QFknDUI/edit#heading=h.46icuqlj57dy
#  AUX2: run each test as its own test script, eg "VLASS_mosaic_stakeholder_test_script.py"
#   not done (is this actually necessary? what is the goal of these scripts such that individual files are needed?)
#
##########################################################################
#
#  Tests
#
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
from datetime import datetime

from casatools import imager
from casatasks import casalog, impbcor, imdev, imhead, imsubimage, imstat, immath
from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
from casatestutils.stakeholder import StkUnitTest

quick_test = False if ('QUICK_TEST' not in os.environ) else os.environ['QUICK_TEST']
quick_test = True if str(quick_test).lower() in ['1', 'true'] else False
use_partial_results = False if ('USE_PARTIAL_RESULTS' not in os.environ) else os.environ['USE_PARTIAL_RESULTS']
use_partial_results = True if str(use_partial_results).lower() in ['1', 'true'] else False
if (quick_test):
    casalog.post("QUICK_TEST env variable found\nRunning tests with reduced image sizes", "INFO")
else:
    casalog.post("QUICK_TEST env variable not found or is false\nRunning tests with full image sizes", "INFO")
casalog.post(f"USE_PARTIAL_RESULTS: {use_partial_results}", "INFO")

##############################################
##############################################
class test_j1302(StkUnitTest):

    def setUp(self):
        super().setUp()
        self.vis = 'J1302-12fields.ms'
        self.phasecenter = '13:03:13.874 -10.51.16.73'
        self.im = imager()

    # Test 1
    def test_j1302_mtmfs(self):
        """ [j1302] test_j1302_mtmfs """
        ######################################################################################
        # Should match values for "Stokes I" in the "Values to be compared"
        ######################################################################################
        # "noncube" to allow give this test a unique prefix, for running with runtest
        # For example: runtest.py -v test_vlass_1v2.py[test_j1302_mtmfs]

        # not part of the jupyter scripts
        tstobj = self # jupyter equivalent: "tstobj = test_j1302()"

        #######################################################
        # %% Set local vars [test_j1302_mtmfs] start @
        #######################################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        data_path_dir = 'J1302/Stakeholder-test-mosaic-data'
        img0 = 'J1302_iter2'
        tstobj.prepData(self.vis, data_path_dir, 'secondmask.mask','QLcatmask.mask')
        imsize = 4000
        spw = ''
        rms = [0.00017975829898762892, 0.0013099727978948515] # tt0, tt1 noise floor as measured from a full-scale image run, Range: [700,800],[3300,1900]

        #####################################################
        # %% Set local vars [test_j1302_mtmfs] end @
        # .......................................

        # not part of the jupyter scripts
        starttime = datetime.now()
        if quick_test:
            imsize = 1000
            tstobj.resize_mask('QLcatmask.mask', 'QLcatmask.mask', [imsize, imsize])
            tstobj.resize_mask('secondmask.mask', 'secondmask.mask', [imsize, imsize])

        # .......................................
        # %% Prepare masks [test_j1302_mtmfs] start @
        ######################################################

        # combine first and 2nd order masks
        if not use_partial_results:
            immath(imagename=['secondmask.mask','QLcatmask.mask'],expr='IM0+IM1',outfile='sum_of_masks.mask')
            tstobj.im.mask(image='sum_of_masks.mask',mask='combined.mask',threshold=0.5)

        ####################################################
        # %% Prepare masks [test_j1302_mtmfs] end @
        # %% Run tclean [test_j1302_mtmfs] start  @
        ####################################################

        records = []
        def run_tclean(vis=tstobj.vis, intent='OBSERVE_TARGET#UNSPECIFIED', uvrange='<12km',
                       niter=None, compare_tclean_pars=None, datacolumn=None,
                       imagename=img0, phasecenter=tstobj.phasecenter, reffreq='3.0GHz',
                       mask='', usemask='user', pbmask=0.0, deconvolver='mtmfs',
                       cell='0.6arcsec', imsize=imsize, spw=spw, gridder='mosaic',
                       uvtaper=[''], restoringbeam=[], mosweight=False, rotatepastep=5.0,
                       smallscalebias=0.4, pblimit=0.1, scales=[0], weighting='briggs',
                       robust=1.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0,
                       savemodel='none', interactive=0, calcres=False, calcpsf=False):
            params = locals()
            params = {k:params[k] for k in filter(lambda x: x not in ['tstobj', 'records'], params.keys())}
            records.append( tstobj.run_tclean(**params) )
            return records[-1]

        script_pars_vals_0 = tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False)
        script_pars_vals_1 = tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)
        script_pars_vals_2 = tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='J1302_iter2', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='modelcolumn', calcres=False, calcpsf=False, parallel=False)
        script_pars_vals_3 = tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)
        script_pars_vals_4 = tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)

        # initialize iter2, no cleaning
        run_tclean( niter=0,     datacolumn='corrected', calcres=True, calcpsf=True,                         compare_tclean_pars=script_pars_vals_0 )

        # # resume iter2 with QL mask
        run_tclean( niter=20000, datacolumn='corrected', mask="QLcatmask.mask", nsigma=3.0, scales=[0,5,12], compare_tclean_pars=script_pars_vals_1 )

        # save model column, doesn't happen here in acutal VLASS pipeline, but makes sure functionality works.
        run_tclean( niter=0,     datacolumn='data',      savemodel='modelcolumn',                            compare_tclean_pars=script_pars_vals_2 )

        # resume iter2 with combined mask, remove old mask first, pass new mask as parameter
        os.system(f"rm -rf {img0}.mask")
        run_tclean( niter=20000, datacolumn='corrected', mask="combined.mask",  nsigma=3.0, scales=[0,5,12], compare_tclean_pars=script_pars_vals_3 )

        # resume iter2 with pbmask, removed old mask first then specify pbmask in resumption of tclean
        os.system(f"rm -rf {img0}.mask")
        run_tclean( niter=20000, datacolumn='corrected', usemask='pb', mask="", pbmask=0.4,   nsigma=4.5, scales=[0,5,12], cycleniter=100, compare_tclean_pars=script_pars_vals_4 )

        ################################################################
        # %% Run tclean [test_j1302_mtmfs] end                @
        # %% Compare Expected Values [test_j1302_mtmfs] start @
        ################################################################

        stats613 = {
            'full': {
                'tt0': 0.3198292,
                'tt1': 0.01994022,
                'alpha': 0.06234646,
                'beam': { 'maj': 3.1565470695495605, 'min': 2.58677792549133, 'pos': 11.282347679138184 },
                'runtime': 2413
            },
            '1/4': {
                'tt0': 0.29676201939582825,
                'tt1': -0.004495048895478249,
                'alpha': -0.015146981924772263,
                'beam': { 'maj': 3.138805627822876, 'min': 2.5604472160339355, 'pos': 10.93341064453125 },
                'runtime': 417
            }
        }
        stats613 = stats613['full'] if (not quick_test) else stats613[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        # N/A for this test

        halfsize   = round(imsize / 4000 * 2000)
        box        = f"{halfsize},{halfsize},{halfsize},{halfsize}"
        tt0stats   = imstat(imagename=img0+'.image.tt0', box=box)
        tt1stats   = imstat(imagename=img0+'.image.tt1', box=box)
        alphastats = imstat(imagename=img0+'.alpha', box=box)
        curr_stats    = np.squeeze(np.array([ tt0stats['max'], tt1stats['max'], alphastats['max']]))
        onaxis_stats  = np.array([            0.3337,         -0.01588,        -0.0476])
        casa613_stats = np.array([            stats613['tt0'], stats613['tt1'], stats613['alpha']])

        # (a) tt0 vs 6.1.3, on-axis
        success0, report0 = tstobj.check_metrics_flux(curr_stats[0], onaxis_stats[0],  valname="Frac Diff tt0 vs. on-axis", rms_or_std=rms[0])
        success1, report1 = tstobj.check_metrics_flux(curr_stats[0], casa613_stats[0], valname="Frac Diff tt0 vs. 6.1.3 image", rms_or_std=rms[0])

        # (b) tt1 vs 6.1.3, on-axis
        success2, report2 = tstobj.check_metrics_flux(curr_stats[1], onaxis_stats[1],  valname="Frac Diff tt1 vs. on-axis", rms_or_std=rms[1])
        success3, report3 = tstobj.check_metrics_flux(curr_stats[1], casa613_stats[1], valname="Frac Diff tt1 vs. 6.1.3 image", rms_or_std=rms[1])

        # (c) alpha images
        success4, report4 = tstobj.check_metrics_alpha(curr_stats[2], onaxis_stats[2],  valname="Abs Diff alpha vs. on-axis", rmss_or_stds=rms)
        success5, report5 = tstobj.check_metrics_alpha(curr_stats[2], casa613_stats[2], valname="Abs Diff alpha vs. 6.1.3 image", rmss_or_stds=rms)

        # (d) beamsize comparison vs 6.1.3
        restbeam          = imhead(img0+'.image.tt0')['restoringbeam']
        beamstats_curr    = np.array([restbeam['major']['value'], restbeam['minor']['value'], restbeam['positionangle']['value']])
        beamstats_613     = np.array([stats613['beam']['maj'],    stats613['beam']['min'],    stats613['beam']['pos']])
        success6, report6 = tstobj.check_fracdiff(beamstats_curr, beamstats_613, valname="Frac Diff Maj, Min, PA vs 6.1.3")

        # (e) Confirm presence of model column in resultant MS
        success7, report7 = tstobj.check_column_exists("MODEL_DATA")

        report  = "".join([report0, report1, report2, report3, report4, report5, report6, report7])
        success = success0 and success1 and success2 and success3 and success4 and success5 and success6 and success7 and tstobj.th.check_final(report)
        casalog.post(f"{report}\nSuccess: {success}", "INFO")

        ##############################################################
        # %% Compare Expected Values [test_j1302_mtmfs] end @
        # %% Generate Images [test_j1302_mtmfs] start       @
        ##############################################################

        if quick_test:
            # self.mom8_creator(image=img0+'.image.tt0', range_list=[0, 0.32])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
            # self.mom8_creator(image=img0+'.image.tt1', range_list=[0, 0.05])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
        else:
            # self.mom8_creator(image=img0+'.image.tt0', range_list=[0, 0.30])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
            # self.mom8_creator(image=img0+'.image.tt1', range_list=[0, 0.01])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")

        ######################################################
        # %% Generate Images [test_j1302_mtmfs] end @
        ######################################################
        # not part of the jupyter scripts

        # save results for future analysis
        np.save(tstobj.id()+'.tt0tt1alpha.npy', curr_stats)
        np.save(tstobj.id()+'.beamstats.npy', beamstats_curr)
        np.save(tstobj.id()+'.tcleanrecs.npy', records)

        # (f) Runtimes not significantly different relative to previous runs
        # don't test this in jupyter notebooks - runtimes differ too much between machines
        if ParallelTaskHelper.isMPIEnabled():
            # runtime with MPI -n 8
            successr, reportr = tstobj.check_runtime(starttime, 5373)
        else:
            successr, reportr = tstobj.check_runtime(starttime, stats613['runtime'])
        report += reportr
        success = success and successr and self.th.check_final(report)

        tstobj.assertTrue(success, msg=report)

    # Test 2
    def test_j1302_awproject(self):
        """ [j1302] test_j1302_awproject """
        ######################################################################################
        # Should match values for "Stokes I" in the "Values to be compared"
        ######################################################################################
        # not part of the jupyter scripts
        tstobj = self # jupyter equivalent: "tstobj = test_j1302()"

        ##################################################
        # %% Set local vars [test_j1302_awproject] start @
        ##################################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        data_path_dir  = 'J1302/Stakeholder-test-awproject-data'
        img0 = 'J1302_iter0d'
        img1 = 'J1302_iter2'
        cache0name, cache1name = "cache0d.cf", "cache2.cf"
        if quick_test:
            tstobj.prepData(self.vis, data_path_dir, f"cfcache_quick1/{img0}.cf", f"cfcache_quick1/{img1}.cf", "QLcatmask.mask", "secondmask.mask")
            os.system(f"mv cfcache_quick1/{img0}.cf {cache0name}")
            os.system(f"mv cfcache_quick1/{img1}.cf {cache1name}")
        else:
            tstobj.prepData(self.vis, data_path_dir, f"cfcache/{img0}.cf", f"cfcache/{img1}.cf", "QLcatmask.mask", "secondmask.mask")
            os.system(f"mv cfcache/{img0}.cf {cache0name}")
            os.system(f"mv cfcache/{img1}.cf {cache1name}")
            # cache0name, cache1name = "", ""
        imsize = 5250
        rms = [0.00025982361923319354, 0.00211483438886223] # tt0, tt1 noise floor as measured from a full-scale image run, Range: [1500,500],[3500,2000]

        ################################################
        # %% Set local vars [test_j1302_awproject] end @
        # .......................................

        # not part of the jupyter scripts
        starttime = datetime.now()
        if quick_test:
            imsize = 1312
            tstobj.resize_mask('QLcatmask.mask', 'QLcatmask.mask', [imsize,imsize])
            tstobj.resize_mask('secondmask.mask', 'secondmask.mask', [imsize,imsize])

        # .......................................
        # %% Prepare masks [test_j1302_mtmfs] start @
        ######################################################

        # combine first and 2nd order masks
        if not use_partial_results:
            immath(imagename=['secondmask.mask','QLcatmask.mask'],expr='IM0+IM1',outfile='sum_of_masks.mask')
            tstobj.im.mask(image='sum_of_masks.mask',mask='combined.mask',threshold=0.5)

        ####################################################
        # %% Prepare masks [test_j1302_mtmfs] end @
        # %% Run tclean [test_j1302_awproject] start       @
        ####################################################

        records = []
        def run_tclean(vis=tstobj.vis, datacolumn='corrected', intent='OBSERVE_TARGET#UNSPECIFIED',
                       imagename=None, niter=None, compare_tclean_pars=None,
                       phasecenter=tstobj.phasecenter, reffreq='3.0GHz', deconvolver='mtmfs',
                       cell='0.6arcsec', imsize=imsize, mask='', usemask='user', pbmask=0.0, cfcache="",
                       wbawp=True, gridder='awproject', uvtaper=[''], uvrange='<12km',
                       restoringbeam=[], wprojplanes=32, mosweight=False, conjbeams=True,
                       usepointing=True, rotatepastep=5.0, pblimit=0.02, scales=[0],
                       weighting='briggs', robust=1.0, nsigma=2.0, cycleniter=5000,
                       savemodel='none', calcpsf=True, calcres=True, cyclefactor=3, interactive=0,
                       smallscalebias=0.4, pointingoffsetsigdev=[300, 30]):
            params = locals()
            params = {k:params[k] for k in filter(lambda x: x not in ['tstobj', 'records'], params.keys())}
            records.append( tstobj.run_tclean(**params) )
            return records[-1]

        def replace_psf(old, new):
            """ Replaces [old] PSF image with [new] image. Clears parallel working directories."""
            for this_tt in ['tt0', 'tt1', 'tt2']:
                shutil.rmtree(old+'.psf.'+this_tt)
                shutil.copytree(new+'.psf.'+this_tt, old+'.psf.'+this_tt)

        script_pars_vals_0 = tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter0d', imsize=5250, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='awproject', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=32, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=False, conjbeams=True, cfcache='', usepointing=True, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[300, 30], pblimit=0.02, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=5000, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=True, parallel=False )
        script_pars_vals_1 = tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=5250, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='awproject', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=32, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=True, cfcache='', usepointing=True, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[300, 30], pblimit=0.02, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=5000, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False )
        script_pars_vals_2 = tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=5250, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='awproject', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=32, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=True, cfcache='', usepointing=True, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[300, 30], pblimit=0.02, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=3000, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False )
        script_pars_vals_3 = tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='J1302_iter2', imsize=5250, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='awproject', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=32, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=True, cfcache='', usepointing=True, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[300, 30], pblimit=0.02, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=5000, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='modelcolumn', calcres=False, calcpsf=False, parallel=False )
        script_pars_vals_4 = tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=5250, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='awproject', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=32, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=True, cfcache='', usepointing=True, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[300, 30], pblimit=0.02, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=3000, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False )
        script_pars_vals_5 = tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2', imsize=5250, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='awproject', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=32, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=True, cfcache='', usepointing=True, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[300, 30], pblimit=0.02, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False )

        # create robust psfs with wbawp=False using corrected column
        run_tclean(imagename=img0, cfcache=cache0name, niter=0, calcres=False, wbawp=False,
                   compare_tclean_pars=script_pars_vals_0)

        # initialize iter2, no cleaning
        run_tclean(imagename=img1, cfcache=cache1name, niter=0, compare_tclean_pars=script_pars_vals_1)

        #  replace iter2 psf with no_WBAP
        replace_psf(img1, img0)

        #  resume iter2 with QL mask
        run_tclean(imagename=img1, cfcache=cache1name, niter=20000, scales=[0, 5, 12], nsigma=3.0, cycleniter=3000,
                   mask="QLcatmask.mask", calcres=False, calcpsf=False, compare_tclean_pars=script_pars_vals_2)

        # save model column, doesn't happen here in acutal VLASS pipeline, but makes sure functionality works.
        run_tclean(imagename=img1, cfcache=cache1name, niter=0, datacolumn='data', calcres=False, calcpsf=False,
                   savemodel='modelcolumn', compare_tclean_pars=script_pars_vals_3)

        # resume iter2 with combined mask
        os.system(f"rm -rf {img1}.mask")
        run_tclean(imagename=img1, cfcache=cache1name, niter=20000, scales=[0, 5, 12], nsigma=3.0, cycleniter=3000,
                   mask="combined.mask", calcres=False, calcpsf=False, compare_tclean_pars=script_pars_vals_4)

        # resume iter2 with pbmask, removed old mask first then specify pbmask in resumption of tclean
        os.system(f"rm -rf {img1}.mask")
        run_tclean(imagename=img1, cfcache=cache1name, niter=20000, scales=[0, 5, 12], nsigma=4.5, cycleniter=500,
                   mask="", calcres=False, calcpsf=False, usemask='pb', pbmask=0.4, compare_tclean_pars=script_pars_vals_5)

        ###########################################################
        # %% Run tclean [test_j1302_awproject] end                @
        # %% Compare Expected Values [test_j1302_awproject] start @
        ###########################################################

        stats613 = {
            'full': {
                'tt0': 0.3174496,
                'tt1': -0.01514572,
                'alpha': -0.04771062,
                'beam': { 'maj': 3.07221413,        'min': 2.49312615,      'pos': 11.04310322 },
                'runtime': 21179
            },
            '1/4': {
                'tt0': 0.21771885454654694,
                'tt1': -0.05718548595905304,
                'alpha': -0.26265749335289,
                'beam': { 'maj': 3.091102123260498, 'min': 2.5045325756073, 'pos': 11.218000411987305 },
                'runtime': 1993 # this is with cfcache
            }
        }
        stats613 = stats613['full'] if (not quick_test) else stats613[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        # N/A: no pbcore, rms, or subim images are created for this test

        halfsize   = round(imsize / 5250 * 2625)
        box        = f"{halfsize},{halfsize},{halfsize},{halfsize}"
        tt0stats=imstat(  imagename=img1+'.image.tt0',box=box)
        tt1stats=imstat(  imagename=img1+'.image.tt1',box=box)
        alphastats=imstat(imagename=img1+'.alpha',    box=box)
        curr_stats=np.squeeze(np.array([tt0stats['max'], tt1stats['max'], alphastats['max']]))
        onaxis_stats=np.array([         0.3337,         -0.01588,        -0.0476])
        casa613_stats=np.array([        stats613['tt0'], stats613['tt1'], stats613['alpha']])

        # (a) tt0 vs 6.1.3, on-axis
        success0, report0 = tstobj.check_metrics_flux(curr_stats[0], onaxis_stats[0],  valname="Frac Diff tt0 vs. on-axis", rms_or_std=rms[0])
        success1, report1 = tstobj.check_metrics_flux(curr_stats[0], casa613_stats[0], valname="Frac Diff tt0 vs. 6.1.3 image", rms_or_std=rms[0])

        # (b) tt1 vs 6.1.3, on-axis
        success2, report2 = tstobj.check_metrics_flux(curr_stats[1], onaxis_stats[1],  valname="Frac Diff tt1 vs. on-axis", rms_or_std=rms[1])
        success3, report3 = tstobj.check_metrics_flux(curr_stats[1], casa613_stats[1], valname="Frac Diff tt1 vs. 6.1.3 image", rms_or_std=rms[1])

        # (c) alpha images
        success4, report4 = tstobj.check_metrics_alpha(curr_stats[2], onaxis_stats[2],  valname="Abs Diff alpha vs. on-axis", rmss_or_stds=rms)
        success5, report5 = tstobj.check_metrics_alpha(curr_stats[2], casa613_stats[2], valname="Abs Diff alpha vs. 6.1.3 image", rmss_or_stds=rms)

        # (d) beamsize comparison vs 6.1.3
        restbeam          = imhead(img1+'.image.tt0')['restoringbeam']
        beamstats_curr    = np.array([restbeam['major']['value'], restbeam['minor']['value'], restbeam['positionangle']['value']])
        beamstats_613     = np.array([stats613['beam']['maj'],    stats613['beam']['min'],    stats613['beam']['pos']])
        success6, report6 = tstobj.check_fracdiff(beamstats_curr, beamstats_613, valname="Frac Diff Maj, Min, PA vs 6.1.3")

        # (e) Confirm presence of model column in resultant MS
        success7, report7 = tstobj.check_column_exists("MODEL_DATA")

        report  = "".join([report0, report1, report2, report3, report4, report5, report6, report7])
        success = success0 and success1 and success2 and success3 and success4 and success5 and success6 and success7 and tstobj.th.check_final(report)
        casalog.post(f"{report}\nSuccess: {success}", "INFO")

        #########################################################
        # %% Compare Expected Values [test_j1302_awproject] end @
        # %% Generate Images [test_j1302_awproject] start       @
        #########################################################

        if quick_test:
            # self.mom8_creator(image=img1+'.image.tt0', range_list=[0, 0.22])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
        else:
            # self.mom8_creator(image=img1+'.image.tt0', range_list=[0, 0.32])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")

        #################################################
        # %% Generate Images [test_j1302_awproject] end @
        #################################################
        # not part of the jupyter scripts

        # save results for future analysis
        np.save(tstobj.id()+'.tt0stats.npy', tt0stats)
        np.save(tstobj.id()+'.tt1stats.npy', tt1stats)
        np.save(tstobj.id()+'.alphastats.npy', alphastats)
        np.save(tstobj.id()+'.beamstats.npy', beamstats_curr)
        np.save(tstobj.id()+'.tcleanrecs.npy', records)

        # (f) Runtimes not significantly different relative to previous runs
        # don't test this in jupyter notebooks - runtimes differ too much between machines
        if ParallelTaskHelper.isMPIEnabled():
            # runtime with MPI -n 8
            successr, reportr = tstobj.check_runtime(starttime, 9594)
        else:
            successr, reportr = tstobj.check_runtime(starttime, stats613['runtime'])
        report += reportr
        success = success and successr and self.th.check_final(report)

        tstobj.assertTrue(success, msg=report)

    # Test 3
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Tclean crashes with mpicasa+mosaic gridder+stokes imaging.")
    def test_j1302_mosaic_cube(self):
        """ [j1302] test_j1302_mosaic_cube """
        ######################################################################################
        # Should match values for "Cube" in the "Values to be compared"
        ######################################################################################
        # not part of the jupyter scripts
        tstobj = self # jupyter equivalent: "tstobj = test_j1302()"

        ####################################################
        # %% Set local vars [test_j1302_mosaic_cube] start @
        ####################################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        data_path_dir  = 'J1302/Stakeholder-test-mosaic-cube-data'
        tstobj.prepData(self.vis, data_path_dir, "QLcatmask.mask", "combined.mask")
        imsize = 4000
        spw_chans = ''
        rms = {'0': 0.0005612289083201638, '1': 0.0004997134396517132, '2': 0.0008543526968930547} # per-spw noise floor as measured from a full-scale image run, Range:[100,100],[3900,1900]

        # reference frequence to use per spectral window (spw)
        refFreqDict  = {
            '0' :  '2.028GHz',
            '1' :  '2.796GHz',
            '2' :  '3.564GHz'
        }

        ##################################################
        # %% Set local vars [test_j1302_mosaic_cube] end @
        # .......................................

        # not part of the jupyter scripts
        starttime = datetime.now()
        if quick_test:
            imsize = 1000
            tstobj.resize_mask("QLcatmask.mask", "QLcatmask.mask", [imsize,imsize])
            tstobj.resize_mask("combined.mask", "combined.mask", [imsize,imsize])

        # .......................................
        # %% Run tclean [test_j1302_mosaic_cube] start   @
        ##################################################

        def iname(image_iter, spw, stokes):
            return 'J1302_'+image_iter+'_'+spw.replace('~','-')+'_'+stokes

        def run_tclean(vis=tstobj.vis, uvrange='<12km', imsize=imsize, intent='OBSERVE_TARGET#UNSPECIFIED',
                       imagename=None, niter=None, compare_tclean_pars=None, spw=None,
                       cell='0.6arcsec', datacolumn='data', phasecenter=tstobj.phasecenter,
                       stokes='I', reffreq='3.0GHz', deconvolver='mtmfs', gridder='mosaic',
                       restoringbeam=[], usemask='user', mask='', pbmask=0.0, mosweight=False,
                       rotatepastep=5.0, pblimit=0.1, scales=[0], nterms=1, weighting='briggs',
                       smallscalebias=0.4, robust=1.0, uvtaper=[''], nsigma=2.0, cycleniter=500,
                       cyclefactor=3, calcpsf=True, calcres=True, interactive=0):
            params = locals()
            params = {k:params[k] for k in filter(lambda x: x not in ['tstobj'], params.keys())}
            return tstobj.run_tclean(**params)

        spws         = [ '0','1','2']
        stokesParams = ['IQUV']
        spwstats     = { '0': {'freq': 2.028, 'IQUV': np.array([0.0,0.0,0.0,0.0]), 'beam': np.array([0.0,0.0,0.0])},
                         '1': {'freq': 2.796, 'IQUV': np.array([0.0,0.0,0.0,0.0]), 'beam': np.array([0.0,0.0,0.0])},
                         '2': {'freq': 3.594, 'IQUV': np.array([0.0,0.0,0.0,0.0]), 'beam': np.array([0.0,0.0,0.0])} }

        script_pars_vals_0 = {
            '0': { 'IQUV': tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False) },
            '1': { 'IQUV': tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False) },
            '2': { 'IQUV': tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False) },
        }
        script_pars_vals_1 = {
            '0': { 'IQUV': tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '1': { 'IQUV': tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '2': { 'IQUV': tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
        }
        script_pars_vals_2 = {
            '0': { 'IQUV': tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '1': { 'IQUV': tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '2': { 'IQUV': tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
        }
        script_pars_vals_3 = {
            '0': { 'IQUV': tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '1': { 'IQUV': tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '2': { 'IQUV': tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1302_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
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
                                   stokes=stokes, reffreq=refFreqDict[spw], compare_tclean_pars=script_pars_vals_0[spw][stokes] )

                # resume iter2 with QL mask
                r[1] = run_tclean( imagename=imagename, datacolumn='corrected', scales=[0,5,12], nsigma=3.0, niter=20000, cycleniter=500,
                                   mask="QLcatmask.mask", calcres=False, calcpsf=False, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw], compare_tclean_pars=script_pars_vals_1[spw][stokes] )

                # resume iter2 with combined mask
                os.system('rm -rf *.workdirectory')
                os.system('rm -rf *iter2*.mask')
                r[2] = run_tclean( imagename=imagename, datacolumn='corrected', scales=[0,5,12], nsigma=3.0, niter=20000, cycleniter=500,
                                   mask="combined.mask", calcres=False, calcpsf=False, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw], compare_tclean_pars=script_pars_vals_2[spw][stokes])

                # os.system('rm -rf iter2*.mask')
                r[3] = run_tclean( imagename=imagename, datacolumn='corrected', scales=[0,5,12], nsigma=4.5, niter=20000, cycleniter=100,
                                   mask="", calcres=False, calcpsf=False, usemask='pb', pbmask=0.4, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw], compare_tclean_pars=script_pars_vals_3[spw][stokes] )

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

        stats613 = {
            'full': {
                'F_nu': 0.3127,
                'alpha': 0.04134,
                'runtime': 7167
            },
            '1/4': {
                'F_nu': 0.2835444715643381,
                'alpha': -0.39275161146262283,
                'runtime': 902
            }
        }
        stats613 = stats613['full'] if (not quick_test) else stats613[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        # N/A: no pbcore, rms, or subim images are created for this test

        # (g) Fit F_nu0 and Alpha from three cube planes and compare: 6.1.3, on-axis
        # compare to alpha (ground truth), and 6.1.3 (fitted for spws 2, 8, 14 from mosaic gridder in CASA 6.1.3)
        alpha = popt[0]
        f_nu0 = 10**popt[1]
        curr_stats        = np.squeeze(np.array([f_nu0,            alpha])) # [flux density, alpha]
        onaxis_stats      = np.array([           0.3337,          -0.0476])
        casa613_stats     = np.array([           stats613['F_nu'], stats613['alpha']])
        success0, report0 = tstobj.check_metrics_flux(curr_stats[0], onaxis_stats[0],  valname="Frac Diff F_nu vs. on-axis", rms_or_std=np.mean(list(rms.values())))
        success1, report1 = tstobj.check_metrics_flux(curr_stats[0], casa613_stats[0], valname="Frac Diff F_nu vs. 6.1.3 image", rms_or_std=np.mean(list(rms.values())))
        success2, report2 = tstobj.check_metrics_alpha_fitted(curr_stats[1], onaxis_stats[1],  valname="Abs Diff alpha vs. on-axis", pcov=pcov)
        success3, report3 = tstobj.check_metrics_alpha_fitted(curr_stats[1], casa613_stats[1], valname="Abs Diff alpha vs. 6.1.3 image", pcov=pcov)

        spwstats_613={
            'full': {
                '0': { 'IQUV':    np.array([ 0.3024486,      -0.00169682,     -0.00040808,     -0.00172231]),
                       'rmsIQUV': np.array([ 0.00058096,      0.00053226,      0.00052702,      0.00055715]),
                       'SNR':     np.array([ 520.60473159,    3.18793882,      0.77430771,      3.09130732]),
                       'beam': np.array([ 4.63033438, 3.7626121, 8.42547607]),
                       'freq': 2.028},
                '1': { 'IQUV':    np.array([ 3.24606597e-01, -1.02521779e-04, -2.74724560e-04, -9.82215162e-04]),
                       'rmsIQUV': np.array([ 0.00060875,      0.00056557,      0.00056437,      0.00058742]),
                       'SNR':     np.array([ 5.33232366e+02,  1.81271795e-01,  4.86780379e-01,  1.67209226e+00]),
                       'beam': np.array([ 3.31582212, 2.81106067, 12.23698997]),
                       'freq': 2.796},
                '2': { 'IQUV':    np.array([ 0.30785725,      0.00112553,     -0.00233958,     -0.00072553]),
                       'rmsIQUV': np.array([ 0.00120004,      0.00096465,      0.00095693,      0.00094171]),
                       'SNR':     np.array([ 256.53982683,    1.16677654,      2.44488894,      0.77044591]),
                       'beam':np.array([ 2.06845379, 1.62075114, 8.15380859]),
                       'freq': 3.564}
            },
            '1/4': {
                '0': { 'IQUV':    np.array([ 0.32101429,     -0.00071505,     -0.00053584,     -0.00094223]),
                       'rmsIQUV': np.array([ 0.00067207,      0.00051443,      0.00050779,      0.00054066]),
                       'SNR':     np.array([ 477.64906346,    1.38996781,      1.05523628,      1.74273892]),
                       'beam': np.array([ 4.70432472, 3.76074767, 9.39245605]),
                       'freq': 2.028 },
                '1': { 'IQUV':    np.array([ 3.11888933e-01,  6.64262800e-04,  6.19260682e-05, -7.18935160e-04]),
                       'rmsIQUV': np.array([ 0.00062024,      0.00046568,      0.00046498,      0.00047676]),
                       'SNR':     np.array([ 5.02853683e+02,  1.42643376e+00,  1.33181314e-01,  1.50796914e+00]),
                       'beam': np.array([  3.35105848,  2.79067874, 12.36252499]),
                       'freq': 2.796 },
                '2': { 'IQUV':    np.array([ 0.25428799,      0.00060686,     -0.00169307,     -0.00073489]),
                       'rmsIQUV': np.array([ 0.00096729,      0.00067782,      0.00066765,      0.00065435]),
                       'SNR':     np.array([ 262.88833666,    0.89530696,      2.53586517,      1.1230845 ]),
                       'beam': np.array([ 2.01840663, 1.6579901, 9.50128174]),
                       'freq': 3.594 },
            }
        }
        spwstats_613 = spwstats_613['full'] if (not quick_test) else spwstats_613[f"1/4"]

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
            # (h) IQUV flux densities of all three spws: 6.1.3
            successN, reportN = tstobj.check_metrics_flux(spwstats[spw]['IQUV'], spwstats_613[spw]['IQUV'],    valname=f"Stokes Comparison (spw {spw}), Frac Diff IQUV vs 6.1.3", rms_or_std=np.mean(list(rms.values())))
            success4.append(successN)
            report4.append(reportN)
            # (i) IQUV flux densities of all three spws: on-axis measurements
            successN, reportN = tstobj.check_metrics_flux(spwstats[spw]['IQUV'], spwstats_onaxis[spw]['IQUV'], valname=f"Stokes Comparison (spw {spw}), Frac Diff IQUV vs on-axis", rms_or_std=np.mean(list(rms.values())))
            success4.append(successN)
            report4.append(reportN)
            # (j) Beam of all three spws:                6.1.3
            successN, reportN = tstobj.check_fracdiff(spwstats[spw]['beam'], spwstats_613[spw]['beam'],        valname=f"Stokes Comparison (spw {spw}), Frac Diff Beam vs 6.1.3")
            success4.append(successN)
            report4.append(reportN)

        report  = "".join([report0, report1, report2, report3, *report4])
        success = success0 and success1 and success2 and success3 and all(success4) and tstobj.th.check_final(report)
        casalog.post(f"{report}\nSuccess: {success}", "INFO")

        ###########################################################
        # %% Compare Expected Values [test_j1302_mosaic_cube] end @
        # %% Generate Images [test_j1302_mosaic_cube] start       @
        ###########################################################

        if quick_test:
            # self.mom8_creator(image=iname('iter2', '0', 'IQUV')+'.image.tt0', range_list=[0, 0.32])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
            # self.mom8_creator(image=iname('iter2', '1', 'IQUV')+'.image.tt0', range_list=[0, 0.32])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
            # self.mom8_creator(image=iname('iter2', '2', 'IQUV')+'.image.tt0', range_list=[0, 0.26])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
        else:
            # self.mom8_creator(image=iname('iter2', '0', 'IQUV')+'.image.tt0', range_list=[0, 0.31])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
            # self.mom8_creator(image=iname('iter2', '1', 'IQUV')+'.image.tt0', range_list=[0, 0.33])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
            # self.mom8_creator(image=iname('iter2', '2', 'IQUV')+'.image.tt0', range_list=[0, 0.31])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")

        ###################################################
        # %% Generate Images [test_j1302_mosaic_cube] end @
        ###################################################
        # not part of the jupyter scripts

        # save results for future analysis
        np.save(tstobj.id()+'.fluxdens_alpha.npy', curr_stats)
        np.save(tstobj.id()+'.freqs.npy', freqs)
        np.save(tstobj.id()+'.fluxes.npy', fluxes)
        np.save(tstobj.id()+'.spwstats.npy', spwstats)
        np.save(tstobj.id()+'.tcleanrecs.npy', records)

        # (f) Runtimes not significantly different relative to previous runs
        # don't test this in jupyter notebooks - runtimes differ too much between machines
        successr, reportr = tstobj.check_runtime(starttime, stats613['runtime'])
        report += reportr
        success = success and successr and self.th.check_final(report)

        tstobj.assertTrue(success, msg=report)

    # Test 4
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Only run in serial, since John Tobin only executed this test in serial (see 01/12/22 comment on CAS-12427).")
    def test_j1302_ql(self):
        """ [j1302] test_j1302_ql """
        ######################################################################################
        # Should match values for "QL" in the "Values to be compared"
        ######################################################################################
        # not part of the jupyter scripts
        tstobj = self # jupyter equivalent: "tstobj = test_j1302()"

        ###########################################
        # %% Set local vars [test_j1302_ql] start @
        ###########################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        data_path_dir  = 'J1302/Stakeholder-test-mosaic-data'
        img0 = 'VLASS1.2.ql.T08t20.J1302.10.2048.v1.I.iter0'
        img1 = 'VLASS1.2.ql.T08t20.J1302.10.2048.v1.I.iter1'
        tstobj.prepData(self.vis, data_path_dir)
        imsize = 7290
        rms = 0.00034846254286391285 # noise floor as measured from a full-scale image run, Range: [3000,3000],[6990,3600]

        #########################################
        # %% Set local vars [test_j1302_ql] end @
        # .......................................

        # not part of the jupyter scripts
        starttime = datetime.now()
        if quick_test:
            imsize = 1822

        # .......................................
        # %% Run tclean [test_j1302_ql] start   @
        #########################################

        records = []
        def run_tclean(vis=tstobj.vis, intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data',
                       imagename=None, niter=None, restoration=None, compare_tclean_pars=None, nsigma=0.0,
                       phasecenter=tstobj.phasecenter, reffreq='3.0GHz', deconvolver='mtmfs',
                       cell='1.0arcsec', imsize=imsize, gridder='mosaic', cycleniter=-1,
                       cyclefactor=1.0, restoringbeam='common', perchanweightdensity=False,
                       mosweight=False, calcres=True, calcpsf=True, scales=[0],
                       weighting='briggs', robust=1.0, interactive=0):
            params = locals()
            params = {k:params[k] for k in filter(lambda x: x not in ['tstobj', 'records'], params.keys())}
            records.append( tstobj.run_tclean(**params) )
            return records[-1]

        script_pars_vals_0 = tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='VLASS1.2.ql.T08t20.J1302.10.2048.v1.I.iter0', imsize=[7290, 7290], cell='1.0arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=False, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=360.0, pointingoffsetsigdev=[], pblimit=0.2, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.0, restoration=False, restoringbeam='common', pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[], niter=0, gain=0.1, threshold='0.0mJy', nsigma=0.0, cycleniter=-1, cyclefactor=1.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=0, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False)
        script_pars_vals_1 = tstobj.get_params_as_dict(vis='J1302-12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='VLASS1.2.ql.T08t20.J1302.10.2048.v1.I.iter1', imsize=[7290, 7290], cell='1.0arcsec', phasecenter='13:03:13.874 -10.51.16.73', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=False, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=360.0, pointingoffsetsigdev=[], pblimit=0.2, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.0, restoration=True, restoringbeam='common', pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=500, cyclefactor=2.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=0, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)
        run_tclean(imagename=img0, niter=0,     restoration=False,                                                                           compare_tclean_pars=script_pars_vals_0)
        if not use_partial_results:
            for ext in ['.weight.tt2', '.weight.tt0', '.psf.tt0', '.residual.tt0', '.weight.tt1', '.sumwt.tt2', '.psf.tt1', '.residual.tt1', '.psf.tt2', '.sumwt.tt1', '.model.tt0', '.pb.tt0', '.model.tt1', '.sumwt.tt0']:
                shutil.copytree(src=img0+ext, dst=img1+ext)
        run_tclean(imagename=img1, niter=20000, restoration=True, nsigma=4.5, cycleniter=500, cyclefactor=2.0, calcres=False, calcpsf=False, compare_tclean_pars=script_pars_vals_1)

        ###########################################
        # %% Run tclean [test_j1302_ql] end       @

        # not part of the jupyter scripts
        if use_partial_results:
            os.system("rm -rf *.pbcor.tt0 *.subim")

        # %% Prepare Images [test_j1302_ql] start @
        ###########################################

        # hifv_pbcor(pipelinemode="automatic")
        for fromext,toext in [('.image.tt0','.image.pbcor.tt0'), ('.residual.tt0','.image.residual.pbcor.tt0')]:
            impbcor(imagename=img1+fromext, pbimage=img1+'.pb.tt0', outfile=img1+toext, mode='divide', cutoff=-1.0, stretch=False)
            tstobj.check_img_exists(img1+toext)

        # hif_makermsimages(pipelinemode="automatic")
        imdev(imagename=img1+'.image.pbcor.tt0',
              outfile=img1+'.image.pbcor.tt0.rms',
              overwrite=True, stretch=False, grid=[10, 10], anchor='ref',
              xlength='60arcsec', ylength='60arcsec', interp='cubic', stattype='xmadm',
              statalg='chauvenet', zscore=-1, maxiter=-1)
        tstobj.check_img_exists(img1+'.image.pbcor.tt0.rms')

        # hif_makecutoutimages(pipelinemode="automatic")
        blc = round(imsize / 7290 * 1785)
        urc = round(imsize / 7290 * 5506)
        for ext in ['.image.tt0', '.residual.tt0', '.image.pbcor.tt0', '.image.pbcor.tt0.rms', '.psf.tt0', '.image.residual.pbcor.tt0', '.pb.tt0']:
            imhead(imagename=img1+ext)
            imsubimage(imagename=img1+ext, outfile=img1+ext+'.subim', box=f"{blc},{blc},{urc},{urc}")
            tstobj.check_img_exists(img1+ext+'.subim')

        ####################################################
        # %% Prepare Images [test_j1302_ql] end            @
        # %% Compare Expected Values [test_j1302_ql] start @
        ####################################################

        stats613 = {
            'full': {
                'F_nu': 0.320879,
                'beam': { 'min': 3.16884375,         'maj': 2.59194756,         'pos': 11.36847878 },
                'runtime': 1805
            },
            '1/4': {
                'F_nu': 0.31872469186782837,
                'beam': { 'min': 3.1795637607574463, 'maj': 2.5885016918182373, 'pos': 11.082575798034668 },
                'runtime': 383
            }
        }
        stats613 = stats613['full'] if (not quick_test) else stats613[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        success0, report0 = tstobj.get_imgs_exist_results()

        # (a) tt0 vs 6.1.3, on-axis
        halfsize          = round(imsize / 7290 * 1860)
        imstat_vals       = imstat(imagename=img1+'.image.pbcor.tt0.subim',box=f"{halfsize},{halfsize},{halfsize},{halfsize}")
        curr_stats        = np.squeeze(np.array([imstat_vals['max']]))
        onaxis_stats      = np.array([           0.3337])
        casa613_stats     = np.array([           stats613['F_nu']])
        success1, report1 = tstobj.check_metrics_flux(curr_stats, onaxis_stats, valname="Frac Diff F_nu vs. on-axis", rms_or_std=rms)
        success2, report2 = tstobj.check_metrics_flux(curr_stats, casa613_stats, valname="Frac Diff F_nu vs. 6.1.3 image", rms_or_std=rms)

        # (b) tt1 vs 6.1.3, on-axis
        # no tt1 images for this test, skip

        # (c) alpha images
        # TODO
        # success3, report3 = ...

        # (d) beamsize comparison vs 6.1.3
        restbeam          = imhead(img1+'.image.pbcor.tt0.subim')['restoringbeam']
        beamstats_curr    = np.array([restbeam['major']['value'], restbeam['minor']['value'], restbeam['positionangle']['value']])
        beamstats_613     = np.array([stats613['beam']['maj'],    stats613['beam']['min'],    stats613['beam']['pos']])
        success4, report4 = tstobj.check_fracdiff(beamstats_curr, beamstats_613, valname="Frac Diff Maj, Min, PA vs 6.1.3")

        report  = "".join([report0, report1, report2, report4])
        success = success1 and success2 and success4 and tstobj.th.check_final(report)
        casalog.post(f"{report}\nSuccess: {success}", "INFO")

        ##################################################
        # %% Compare Expected Values [test_j1302_ql] end @
        # %% Generate Images [test_j1302_ql] start       @
        ##################################################

        if quick_test:
            # self.mom8_creator(image=img1+'.image.pbcor.tt0.subim', range_list=[0, 0.32])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
        else:
            # self.mom8_creator(image=img1+'.image.pbcor.tt0.subim', range_list=[0, 0.32])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")

        ##########################################
        # %% Generate Images [test_j1302_ql] end @
        ##########################################
        # not part of the jupyter scripts

        # save results for future analysis
        np.save(tstobj.id()+'.sourceflux.npy', curr_stats)
        np.save(tstobj.id()+'.beamstats.npy', beamstats_curr)
        np.save(tstobj.id()+'.tcleanrecs.npy', records)

        # (f) Runtimes not significantly different relative to previous runs
        # don't test this in jupyter notebooks - runtimes differ too much between machines
        successr, reportr = tstobj.check_runtime(starttime, stats613['runtime'])
        report += reportr
        success = success and successr and self.th.check_final(report)

        tstobj.assertTrue(success, msg=report)


##############################################
##############################################
class test_j1927(StkUnitTest):

    def setUp(self):
        super().setUp()
        self.vis = 'J1927_12fields.ms'
        self.phasecenter = '19:27:30.443 +61.17.32.898'
        self.im = imager()

    # Test 5
    def test_j1927_mtmfs(self):
        """ [j1927] test_j1927_mtmfs """
        ######################################################################################
        # Should match values for "Stokes I" in the "Values to be compared"
        ######################################################################################
        # "noncube" to allow give this test a unique prefix, for running with runtest
        # For example: runtest.py -v test_vlass_1v2.py[test_j1927_mtmfs]

        # not part of the jupyter scripts
        tstobj = self # jupyter equivalent: "tstobj = test_j1927()"

        #######################################################
        # %% Set local vars [test_j1927_mtmfs] start @
        #######################################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        data_path_dir  = 'J1927/J1927-stakeholdertest-mosaic-data'
        img0 = 'J1927_iter2'
        tstobj.prepData(self.vis, data_path_dir, "secondmask.mask", "QLcatmask.mask")
        imsize = 4000
        spw = ''
        rms = [0.0001483304420688553, 0.0007968044018725578] # tt0, tt1 noise floor as measured from a full-scale image run, Range: [500,500],[3400,1900]

        #####################################################
        # %% Set local vars [test_j1927_mtmfs] end @
        # ...................................................

        # not part of the jupyter scripts
        starttime = datetime.now()
        if quick_test:
            imsize = 1000
            tstobj.resize_mask('QLcatmask.mask', 'QLcatmask.mask', [imsize,imsize])
            tstobj.resize_mask('secondmask.mask', 'secondmask.mask', [imsize,imsize])

        # ....................................................
        # %% Prepare masks [test_j1927_mtmfs] start @
        ######################################################

        # combine first and 2nd order masks
        if not use_partial_results:
            immath(imagename=['secondmask.mask','QLcatmask.mask'],expr='IM0+IM1',outfile='sum_of_masks.mask')
            tstobj.im.mask(image='sum_of_masks.mask',mask='combined.mask',threshold=0.5)

        ####################################################
        # %% Prepare masks [test_j1927_mtmfs] end @
        # %% Run tclean [test_j1927_mtmfs] start  @
        ####################################################

        records = []
        def run_tclean(vis=tstobj.vis, intent='OBSERVE_TARGET#UNSPECIFIED', uvrange='<12km',
                       niter=None, compare_tclean_pars=None, datacolumn=None,
                       imagename=img0, phasecenter=tstobj.phasecenter, reffreq='3.0GHz',
                       mask='', savemodel='none', usemask='user', pbmask=0.0,
                       deconvolver='mtmfs', cell='0.6arcsec', imsize=imsize, spw=spw,
                       gridder='mosaic', uvtaper=[''], restoringbeam=[], mosweight=False,
                       rotatepastep=5.0, smallscalebias=0.4, pblimit=0.1, scales=[0],
                       weighting='briggs', robust=1.0, nsigma=2.0, cycleniter=500,
                       cyclefactor=3.0, interactive=0, calcres=False, calcpsf=False):
            params = locals()
            params = {k:params[k] for k in filter(lambda x: x not in ['tstobj', 'records'], params.keys())}
            records.append( tstobj.run_tclean(**params) )
            return records[-1]

        # These are the parameter values from running tclean with John's scripts on the CAS-12427 ticket
        script_pars_vals_0 = tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False)
        script_pars_vals_1 = tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)
        script_pars_vals_2 = tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='J1927_iter2', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='modelcolumn', calcres=False, calcpsf=False, parallel=False)
        script_pars_vals_3 = tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)
        script_pars_vals_4 = tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=2, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)

        # initialize iter2, no cleaning
        run_tclean( niter=0,     datacolumn='corrected', calcres=True, calcpsf=True,                         compare_tclean_pars=script_pars_vals_0 )

        # # resume iter2 with QL mask
        run_tclean( niter=20000, datacolumn='corrected', mask="QLcatmask.mask", nsigma=3.0, scales=[0,5,12], compare_tclean_pars=script_pars_vals_1 )

        # save model column, doesn't happen here in acutal VLASS pipeline, but makes sure functionality works.
        run_tclean( niter=0,     datacolumn='data',      savemodel='modelcolumn',                            compare_tclean_pars=script_pars_vals_2 )

        # resume iter2 with combined mask, remove old mask first, pass new mask as parameter
        os.system(f"rm -rf {img0}.mask")
        run_tclean( niter=20000, datacolumn='corrected', mask="combined.mask",  nsigma=3.0, scales=[0,5,12], compare_tclean_pars=script_pars_vals_3 )

        # resume iter2 with pbmask, removed old mask first then specify pbmask in resumption of tclean
        os.system(f"rm -rf {img0}.mask")
        run_tclean( niter=20000, datacolumn='corrected', usemask='pb', mask="", pbmask=0.4,   nsigma=4.5, scales=[0,5,12], cycleniter=100, compare_tclean_pars=script_pars_vals_4 )

        ################################################################
        # %% Run tclean [test_j1927_mtmfs] end                @
        # %% Compare Expected Values [test_j1927_mtmfs] start @
        ################################################################

        stats613 = {
            'full': {
                'tt0': 0.8887,
                'tt1': 0.4148,
                'alpha': 0.4668,
                'beam': { 'min': 2.46035814,         'maj': 2.05488539,        'pos': -23.22724915 },
                'runtime': 4415
            },
            '1/4': {
                'tt0': 0.878496527671814,
                'tt1': 0.32516300678253174,
                'alpha': 0.3701357841491699,
                'beam': { 'min': 2.3685483932495117, 'maj': 1.971137285232544, 'pos': -24.036651611328125 },
                'runtime': 436
            }
        }
        stats613 = stats613['full'] if (not quick_test) else stats613[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        # N/A for this test

        halfsize          = round(imsize / 4000 * 2000)
        box               = f"{halfsize},{halfsize},{halfsize},{halfsize}"
        tt0stats   = imstat(imagename=img0+'.image.tt0', box=box)
        tt1stats   = imstat(imagename=img0+'.image.tt1', box=box)
        alphastats = imstat(imagename=img0+'.alpha',     box=box)
        curr_stats    = np.squeeze(np.array([ tt0stats['max'], tt1stats['max'], alphastats['max']]))
        onaxis_stats  = np.array([            0.9509,          0.3601,          0.3796])
        casa613_stats = np.array([            stats613['tt0'], stats613['tt1'], stats613['alpha']])

        # (a) tt0 vs 6.1.3, on-axis
        success0, report0 = tstobj.check_metrics_flux(curr_stats[0], onaxis_stats[0],  valname="Frac Diff tt0 vs. on-axis", rms_or_std=rms[0])
        success1, report1 = tstobj.check_metrics_flux(curr_stats[0], casa613_stats[0], valname="Frac Diff tt0 vs. 6.1.3 image", rms_or_std=rms[0])

        # (b) tt1 vs 6.1.3, on-axis
        success2, report2 = tstobj.check_metrics_flux(curr_stats[1], onaxis_stats[1],  valname="Frac Diff tt1 vs. on-axis", rms_or_std=rms[1])
        success3, report3 = tstobj.check_metrics_flux(curr_stats[1], casa613_stats[1], valname="Frac Diff tt1 vs. 6.1.3 image", rms_or_std=rms[1])

        # (c) alpha images
        success4, report4 = tstobj.check_metrics_alpha(curr_stats[2], onaxis_stats[2],  valname="Abs Diff alpha vs. on-axis", rmss_or_stds=rms)
        success5, report5 = tstobj.check_metrics_alpha(curr_stats[2], casa613_stats[2], valname="Abs Diff alpha vs. 6.1.3 image", rmss_or_stds=rms)

        # (d) beamsize comparison vs 6.1.3
        restbeam          = imhead(img0+'.image.tt0')['restoringbeam']
        beamstats_curr    = np.array([restbeam['major']['value'], restbeam['minor']['value'], restbeam['positionangle']['value']])
        beamstats_613     = np.array([stats613['beam']['maj'],    stats613['beam']['min'],    stats613['beam']['pos']])
        success6, report6 = tstobj.check_fracdiff(beamstats_curr, beamstats_613, valname="Frac Diff Maj, Min, PA vs 6.1.3")

        # (e) Confirm presence of model column in resultant MS
        success7, report7 = tstobj.check_column_exists("MODEL_DATA")

        report  = "".join([report0, report1, report2, report3, report4, report5, report6, report7])
        success = success0 and success1 and success2 and success3 and success4 and success5 and success6 and success7 and tstobj.th.check_final(report)
        casalog.post(f"{report}\nSuccess: {success}", "INFO")

        ##############################################################
        # %% Compare Expected Values [test_j1927_mtmfs] end @
        # %% Generate Images [test_j1927_mtmfs] start       @
        ##############################################################

        if quick_test:
            # self.mom8_creator(image=img0+'.image.tt0', range_list=[0, 0.88])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
            # self.mom8_creator(image=img0+'.image.tt1', range_list=[0, 0.33])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
        else:
            # self.mom8_creator(image=img0+'.image.tt0', range_list=[0, 0.89])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
            # self.mom8_creator(image=img0+'.image.tt1', range_list=[0, 0.42])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")

        ######################################################
        # %% Generate Images [test_j1927_mtmfs] end @
        ######################################################
        # not part of the jupyter scripts

        # save results for future analysis
        np.save(tstobj.id()+'.tt0tt1alpha.npy', curr_stats)
        np.save(tstobj.id()+'.beamstats.npy', beamstats_curr)
        np.save(tstobj.id()+'.tcleanrecs.npy', records)

        # (f) Runtimes not significantly different relative to previous runs
        # don't test this in jupyter notebooks - runtimes differ too much between machines
        if ParallelTaskHelper.isMPIEnabled():
            # runtime with MPI -n 8
            successr, reportr = tstobj.check_runtime(starttime, 7129)
        else:
            successr, reportr = tstobj.check_runtime(starttime, stats613['runtime'])
        report += reportr
        success = success and successr and self.th.check_final(report)

        tstobj.assertTrue(success, msg=report)

    # N/A not implemented
    # def test_j1927_awproject(self):
    #     pass

    # Test 6
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Skip test. Tclean crashes with mpicasa+mosaic gridder+stokes imaging.")
    def test_j1927_mosaic_cube(self):
        """ [j1927] test_j1927_mosaic_cube """
        ######################################################################################
        # Should match values for "Cube" in the "Values to be compared"
        ######################################################################################
        # not part of the jupyter scripts
        tstobj = self # jupyter equivalent: "tstobj = test_j1927()"

        ####################################################
        # %% Set local vars [test_j1927_mosaic_cube] start @
        ####################################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        data_path_dir  = 'J1927/J1927-stakeholdertest-mosaic-cube-data'
        tstobj.prepData(self.vis, data_path_dir, "QLcatmask.mask", "combined.mask")
        # rundir = "/users/bbean/dev/CAS-12427/src/casalith/build-casalith/work/linux/test_vlass_j1927_cube_unittest"
        # os.system(f"mv {rundir}/run_results/VLASS* {rundir}/nosedir/test_vlass_1v2/")
        imsize = 4000
        spw_chans = ''
        rms = {'0': 0.0004637447465635465, '1': 0.0005062325702676312, '2': 0.0004403419438565781} # per-spw noise floor as measured from a full-scale image run, Range:[100,100],[3900,1900]

        # reference frequence to use per spectral window (spw)
        refFreqDict  = {
            '0' :  '2.028GHz',
            '1' :  '2.796GHz',
            '2' :  '3.564GHz'
        }

        ##################################################
        # %% Set local vars [test_j1927_mosaic_cube] end @
        # .......................................

        # not part of the jupyter scripts
        starttime = datetime.now()
        if quick_test:
            imsize = 1000
            tstobj.resize_mask("QLcatmask.mask", "QLcatmask.mask", [imsize,imsize])
            tstobj.resize_mask("combined.mask", "combined.mask", [imsize,imsize])

        # .......................................
        # %% Run tclean [test_j1927_mosaic_cube] start   @
        ##################################################

        def iname(image_iter, spw, stokes):
            return 'J1927_'+image_iter+'_'+spw.replace('~','-')+'_'+stokes

        def run_tclean(vis=tstobj.vis, uvrange='<12km', imsize=imsize, intent='OBSERVE_TARGET#UNSPECIFIED',
                       imagename=None, niter=None, compare_tclean_pars=None, spw=None,
                       cell='0.6arcsec', datacolumn='data', phasecenter=tstobj.phasecenter,
                       stokes='I', mask='', usemask='user', pbmask=0.0, reffreq='3.0GHz',
                       deconvolver='mtmfs', gridder='mosaic', restoringbeam=[], mosweight=False,
                       rotatepastep=5.0, pblimit=0.1, scales=[0], nterms=1, weighting='briggs',
                       smallscalebias=0.4, robust=1.0, uvtaper=[''], nsigma=2.0, calcres=True,
                       calcpsf=True, cycleniter=500, cyclefactor=3, interactive=0):
            params = locals()
            params = {k:params[k] for k in filter(lambda x: x not in ['tstobj'], params.keys())}
            tstobj.run_tclean(**params)

        spws         = [ '0','1','2']
        stokesParams = ['IQUV']
        spwstats     = { '0': {'freq': 2.028, 'IQUV': np.array([0.0,0.0,0.0,0.0]), 'beam': np.array([0.0,0.0,0.0])},
                         '1': {'freq': 2.796, 'IQUV': np.array([0.0,0.0,0.0,0.0]), 'beam': np.array([0.0,0.0,0.0])},
                         '2': {'freq': 3.594, 'IQUV': np.array([0.0,0.0,0.0,0.0]), 'beam': np.array([0.0,0.0,0.0])} }

        script_pars_vals_0 = {
            '0': { 'IQUV': tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False) },
            '1': { 'IQUV': tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False) },
            '2': { 'IQUV': tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=0, gain=0.1, threshold=0.0, nsigma=2.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False) },
        }
        script_pars_vals_1 = {
            '0': { 'IQUV': tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '1': { 'IQUV': tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '2': { 'IQUV': tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='QLcatmask.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
        }
        script_pars_vals_2 = {
            '0': { 'IQUV': tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '1': { 'IQUV': tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '2': { 'IQUV': tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=3.0, cycleniter=500, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='user', mask='combined.mask', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
        }
        script_pars_vals_3 = {
            '0': { 'IQUV': tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='2', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_2_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.028GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '1': { 'IQUV': tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='8', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_8_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='2.796GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
            '2': { 'IQUV': tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='14', timerange='', uvrange='<12km', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='corrected', imagename='J1927_iter2_14_IQUV', imsize=4000, cell='0.6arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='IQUV', projection='SIN', startmodel='', specmode='mfs', reffreq='3.564GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=5.0, pointingoffsetsigdev=[], pblimit=0.1, normtype='flatnoise', deconvolver='mtmfs', scales=[0, 5, 12], nterms=1, smallscalebias=0.4, restoration=True, restoringbeam=[], pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[''], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=100, cyclefactor=3.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=False, usemask='pb', mask='', pbmask=0.4, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False) },
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
                                   stokes=stokes, reffreq=refFreqDict[spw], compare_tclean_pars=script_pars_vals_0[spw][stokes] )

                # # resume iter2 with QL mask
                r[1] = run_tclean( imagename=imagename, datacolumn='corrected', scales=[0,5,12], nsigma=3.0, niter=20000, cycleniter=500,   
                                   mask="QLcatmask.mask", calcres=False, calcpsf=False, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw], compare_tclean_pars=script_pars_vals_1[spw][stokes] )

                # resume iter2 with combined mask
                os.system('rm -rf *.workdirectory')
                os.system('rm -rf *iter2*.mask')
                r[2] = run_tclean( imagename=imagename, datacolumn='corrected', scales=[0,5,12], nsigma=3.0, niter=20000, cycleniter=500,
                                   mask="combined.mask", calcres=False, calcpsf=False, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw], compare_tclean_pars=script_pars_vals_2[spw][stokes])

                # os.system('rm -rf iter2*.mask')
                r[3] = run_tclean( imagename=imagename, datacolumn='corrected', scales=[0,5,12], nsigma=4.5, niter=20000, cycleniter=100,
                                   mask="", calcres=False, calcpsf=False, usemask='pb', pbmask=0.4, spw=spw_str,
                                   stokes=stokes, reffreq=refFreqDict[spw], compare_tclean_pars=script_pars_vals_3[spw][stokes] )

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

        stats613 = {
            'full': {
                'F_nu': 0.90649462,
                'alpha': 0.4127,
                'runtime': 15329
            },
            '1/4': {
                'F_nu': 0.8809638356491677,
                'alpha': 0.2955769626453426,
                'runtime': 1022
            }
        }
        stats613 = stats613['full'] if (not quick_test) else stats613[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        # N/A: no pbcore, rms, or subim images are created for this test

        # (g) Fit F_nu0 and Alpha from three cube planes and compare: 6.1.3, on-axis
        # compare to alpha (ground truth), and 6.1.3 (fitted for spws 2, 8, 14 from mosaic gridder in CASA 6.1.3)
        alpha = popt[0]
        f_nu0 = 10**popt[1]
        curr_stats        = np.squeeze(np.array([f_nu0,            alpha])) # [flux density, alpha]
        onaxis_stats      = np.array([           0.9509,           0.3601])
        casa613_stats     = np.array([           stats613['F_nu'], stats613['alpha']])
        success0, report0 = tstobj.check_metrics_flux(curr_stats[0], onaxis_stats[0],  valname="Frac Diff F_nu vs. on-axis", rms_or_std=np.mean(list(rms.values())))
        success1, report1 = tstobj.check_metrics_flux(curr_stats[0], casa613_stats[0], valname="Frac Diff F_nu vs. 6.1.3 image", rms_or_std=np.mean(list(rms.values())))
        success2, report2 = tstobj.check_metrics_alpha_fitted(curr_stats[1], onaxis_stats[1],  valname="Abs Diff alpha vs. on-axis", pcov=pcov)
        success3, report3 = tstobj.check_metrics_alpha_fitted(curr_stats[1], casa613_stats[1], valname="Abs Diff alpha vs. 6.1.3 image", pcov=pcov)

        spwstats_613= {
            'full': {
                '0': {'freq': 2.028,
                      'IQUV': np.array([ 0.75571942,  0.00576184,  0.00080162, -0.00630026]),
                      'beam': np.array([  3.65375686, 3.0324676,  -22.7224865 ])},
                '1': {'freq': 2.796,
                      'IQUV': np.array([ 0.86450773,  0.00522543, -0.00604387, -0.00717182]),
                      'beam': np.array([  2.66741753, 2.23318672, -26.65664101])},
                '2': {'freq': 3.594,
                      'IQUV': np.array([ 0.95684659, 0.00543807, -0.00958768, -0.00369545]),
                      'beam': np.array([ 2.0738709,  1.70000803, -27.73123169])}
            },
            '1/4': {
                '0': { 'freq': 2.028,
                       'IQUV': np.array([ 0.7791100144386292, 0.006753005087375641,  0.00043800054118037224, 0.005206027068197727]),
                       'beam': np.array([ 3.659905195236206,  3.018415927886963,     -23.957826614379883]) },
                '1': { 'freq': 2.796,
                       'IQUV': np.array([ 0.876945972442627,  0.0049004945904016495, -0.005897939205169678,  0.008470434695482254]),
                       'beam': np.array([ 2.6499712467193604, 2.2087230682373047,    -26.021034240722656]) },
                '2': { 'freq': 3.594,
                       'IQUV': np.array([ 0.9208499193191528, 0.006781525909900665,  -0.012631140649318695,  0.012594047002494335]),
                       'beam': np.array([ 2.068054437637329,  1.6841689348220825,    -27.218856811523438]) }
            }
        }

        spwstats_613 = spwstats_613['full'] if (not quick_test) else spwstats_613[f"1/4"]

        success4 = []
        report4 = []
        for spw in spws:
            # (h) IQUV flux densities of all three spws:              6.1.3
            successN, reportN = tstobj.check_metrics_flux(spwstats[spw]['IQUV'], spwstats_613[spw]['IQUV'], valname=f"Stokes Comparison (spw {spw}), Frac Diff IQUV vs 6.1.3", rms_or_std=np.mean(list(rms.values())))
            success4.append(successN)
            report4.append(reportN)
            # (i) IQUV flux densities of all three spws:              on-axis measurements
            # N/A: no no-axis measurements available in VLASS_mosaic_cube_stakeholder_test_script.py
            # (j) Beam of all three spws:                             6.1.3
            successN, reportN = tstobj.check_fracdiff(spwstats[spw]['beam'], spwstats_613[spw]['beam'],     valname=f"Stokes Comparison (spw {spw}), Frac Diff Beam vs 6.1.3")
            success4.append(successN)
            report4.append(reportN)

        report  = "".join([report0, report1, report2, report3, *report4])
        success = success0 and success1 and success2 and success3 and all(success4) and tstobj.th.check_final(report)
        casalog.post(f"{report}\nSuccess: {success}", "INFO")

        ##############################################################
        # %% Compare Expected Values [test_j1927_mosaic_cube] end @
        # %% Generate Images [test_j1927_mosaic_cube] start       @
        ##############################################################

        if quick_test:
            # self.mom8_creator(image=iname('iter2', '0', 'IQUV')+'.image.tt0', range_list=[0, 0.78])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
            # self.mom8_creator(image=iname('iter2', '1', 'IQUV')+'.image.tt0', range_list=[0, 0.88])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
            # self.mom8_creator(image=iname('iter2', '2', 'IQUV')+'.image.tt0', range_list=[0, 0.92])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
        else:
            # self.mom8_creator(image=iname('iter2', '0', 'IQUV')+'.image.tt0', range_list=[0, 0.76])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
            # self.mom8_creator(image=iname('iter2', '1', 'IQUV')+'.image.tt0', range_list=[0, 0.87])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
            # self.mom8_creator(image=iname('iter2', '2', 'IQUV')+'.image.tt0', range_list=[0, 0.96])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")

        ######################################################
        # %% Generate Images [test_j1927_mosaic_cube] end @
        ######################################################
        # not part of the jupyter scripts

        # save results for future analysis
        np.save(tstobj.id()+'.fluxdens_alpha.npy', curr_stats)
        np.save(tstobj.id()+'.freqs.npy', freqs)
        np.save(tstobj.id()+'.fluxes.npy', fluxes)
        np.save(tstobj.id()+'.spwstats.npy', spwstats)
        np.save(tstobj.id()+'.tcleanrecs.npy', records)

        # (f) Runtimes not significantly different relative to previous runs
        # don't test this in jupyter notebooks - runtimes differ too much between machines
        successr, reportr = tstobj.check_runtime(starttime, stats613['runtime'])
        report += reportr
        success = success and successr and self.th.check_final(report)

        tstobj.assertTrue(success, msg=report)

    # Test 7
    # @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "Only run in serial, since John Tobin only executed this test in serial (see 01/12/22 comment on CAS-12427).")
    def test_j1927_ql(self):
        """ [j1927] test_j1927_ql """
        ######################################################################################
        # Should match values for "QL" in the "Values to be compared"
        ######################################################################################
        # not part of the jupyter scripts
        tstobj = self # jupyter equivalent: "tstobj = test_j1927()"

        ###########################################
        # %% Set local vars [test_j1927_ql] start @
        ###########################################

        #previous steps in the pipeline would have created mask files from catalogs and images that were created as an
        #intermediate pipeline step.
        data_path_dir  = 'J1927/J1927-stakeholdertest-mosaic-data'
        img0 = 'VLASS1.2.ql.T26t15.J1927.10.2048.v1.I.iter0'
        img1 = 'VLASS1.2.ql.T26t15.J1927.10.2048.v1.I.iter1'
        tstobj.prepData(self.vis, data_path_dir)
        # rundir = "/users/bbean/dev/CAS-12427/src/casalith/build-casalith/work/linux/test_vlass_j1927_QL_unittest"
        # os.system(f"mv {rundir}/run_results/VLASS* {rundir}/nosedir/test_vlass_1v2/")
        imsize = 7290
        rms = 0.00023228885125825126 # noise floor as measured from a full-scale image run, Range: [2800,2900],[4300,3600]

        #########################################
        # %% Set local vars [test_j1927_ql] end @
        # .......................................

        # not part of the jupyter scripts
        starttime = datetime.now()
        if quick_test:
            imsize = 1822

        # .......................................
        # %% Run tclean [test_j1927_ql] start   @
        #########################################

        records = []
        def run_tclean(vis=tstobj.vis, intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data',
                       imagename=None, niter=None, restoration=None, compare_tclean_pars=None,
                       phasecenter=tstobj.phasecenter, reffreq='3.0GHz', deconvolver='mtmfs', nsigma=0.0, cycleniter=-1, cyclefactor=1.0,
                       cell='1.0arcsec', imsize=imsize, gridder='mosaic', restoringbeam='common',
                       perchanweightdensity=False, mosweight=False, scales=[0], calcres=True, calcpsf=True,
                       weighting='briggs', robust=1.0, interactive=0):
            params = locals()
            params = {k:params[k] for k in filter(lambda x: x not in ['tstobj', 'records'], params.keys())}
            records.append( tstobj.run_tclean(**params) )
            return records[-1]

        script_pars_vals_0 = tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='VLASS1.2.ql.T26t15.J1927.10.2048.v1.I.iter0', imsize=[7290, 7290], cell='1.0arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=False, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=360.0, pointingoffsetsigdev=[], pblimit=0.2, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.0, restoration=False, restoringbeam='common', pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[], niter=0, gain=0.1, threshold='0.0mJy', nsigma=0.0, cycleniter=-1, cyclefactor=1.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=0, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False)
        script_pars_vals_1 = tstobj.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='VLASS1.2.ql.T26t15.J1927.10.2048.v1.I.iter1', imsize=[7290, 7290], cell='1.0arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=False, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=360.0, pointingoffsetsigdev=[], pblimit=0.2, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.0, restoration=True, restoringbeam='common', pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[], niter=20000, gain=0.1, threshold=0.0, nsigma=4.5, cycleniter=500, cyclefactor=2.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=0, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=False, calcpsf=False, parallel=False)
        run_tclean(imagename=img0, niter=0,     restoration=False,                                                                           compare_tclean_pars=script_pars_vals_0)
        if not use_partial_results:
            for ext in ['.weight.tt2', '.weight.tt0', '.psf.tt0', '.residual.tt0', '.weight.tt1', '.sumwt.tt2', '.psf.tt1', '.residual.tt1', '.psf.tt2', '.sumwt.tt1', '.model.tt0', '.pb.tt0', '.model.tt1', '.sumwt.tt0']:
                shutil.copytree(src=img0+ext, dst=img1+ext)
        run_tclean(imagename=img1, niter=20000, restoration=True, nsigma=4.5, cycleniter=500, cyclefactor=2.0, calcres=False, calcpsf=False, compare_tclean_pars=script_pars_vals_1)

        ###########################################
        # %% Run tclean [test_j1927_ql] end       @

        # not part of the jupyter scripts
        if use_partial_results:
            os.system("rm -rf *.pbcor.tt0* *.subim")

        # %% Prepare Images [test_j1927_ql] start @
        ###########################################

        # hifv_pbcor(pipelinemode="automatic")
        for fromext,toext in [('.image.tt0','.image.pbcor.tt0'), ('.residual.tt0','.image.residual.pbcor.tt0')]:
            impbcor(imagename=img1+fromext, pbimage=img1+'.pb.tt0', outfile=img1+toext, mode='divide', cutoff=-1.0, stretch=False)
            tstobj.check_img_exists(img1+toext)

        # hif_makermsimages(pipelinemode="automatic")
        imdev(imagename=img1+'.image.pbcor.tt0',
              outfile=img1+'.image.pbcor.tt0.rms',
              overwrite=True, stretch=False, grid=[10, 10], anchor='ref',
              xlength='60arcsec', ylength='60arcsec', interp='cubic', stattype='xmadm',
              statalg='chauvenet', zscore=-1, maxiter=-1)
        tstobj.check_img_exists(img1+'.image.pbcor.tt0.rms')

        # hif_makecutoutimages(pipelinemode="automatic")
        blc = round(imsize / 7290 * 1785)
        urc = round(imsize / 7290 * 5506)
        for ext in ['.image.tt0', '.residual.tt0', '.image.pbcor.tt0', '.image.pbcor.tt0.rms', '.psf.tt0', '.image.residual.pbcor.tt0', '.pb.tt0']:
            imhead(imagename=img1+ext)
            imsubimage(imagename=img1+ext, outfile=img1+ext+'.subim', box=f"{blc},{blc},{urc},{urc}")
            tstobj.check_img_exists(img1+ext+'.subim')

        ####################################################
        # %% Prepare Images [test_j1927_ql] end            @
        # %% Compare Expected Values [test_j1927_ql] start @
        ####################################################

        stats613 = {
            'full': {
                'F_nu': [0.90649462],
                'beam': { 'min': 2.5034000873565674, 'maj': 2.0568439960479736, 'pos': -27.06390953063965 },
                'runtime': 4277
            },
            '1/4': {
                'F_nu': 0.8989461064338684,
                'beam': { 'min': 2.487116813659668,  'maj': 2.04219126701355,   'pos': -27.007415771484375 },
                'runtime': 826
            }
        }
        stats613 = stats613['full'] if (not quick_test) else stats613[f"1/4"]

        # (l) Ensure intermediate products exist, pbcor images, RMS image (made by imdev), and cutouts (.subim) from imsubimage
        success0, report0 = tstobj.get_imgs_exist_results()

        # (a) tt0 vs 6.1.3, on-axis
        halfsize          = round(imsize / 7290 * 1860)
        box               = f"{halfsize},{halfsize},{halfsize},{halfsize}"
        imstat_vals       = imstat(imagename=img1+'.image.pbcor.tt0.subim', box=box)
        curr_stats        = np.squeeze(np.array([imstat_vals['max']]))
        onaxis_stats      = np.array([           0.9509])
        casa613_stats     = np.array([           stats613['F_nu']])
        success1, report1 = tstobj.check_metrics_flux(curr_stats, onaxis_stats,  valname="Frac Diff F_nu vs. on-axis", rms_or_std=rms)
        success2, report2 = tstobj.check_metrics_flux(curr_stats, casa613_stats, valname="Frac Diff F_nu vs. 6.1.3 image", rms_or_std=rms)

        # (b) tt1 vs 6.1.3, on-axis
        # no tt1 images for this test, skip

        # (c) alpha images
        # TODO
        # success3, report3 = ...

        # (d) beamsize comparison vs 6.1.3
        restbeam          = imhead(img1+'.image.pbcor.tt0.subim')['restoringbeam']
        beamstats_curr    = np.array([restbeam['major']['value'], restbeam['minor']['value'], restbeam['positionangle']['value']])
        beamstats_613     = np.array([stats613['beam']['maj'],    stats613['beam']['min'],    stats613['beam']['pos']])
        success4, report4 = tstobj.check_fracdiff(beamstats_curr, beamstats_613, valname="Frac Diff Maj, Min, PA vs 6.1.3")

        report  = "".join([report0, report1, report2, report4])
        success = success1 and success2 and success4 and tstobj.th.check_final(report)
        casalog.post(f"{report}\nSuccess: {success}", "INFO")

        ##################################################
        # %% Compare Expected Values [test_j1927_ql] end @
        # %% Generate Images [test_j1927_ql] start       @
        ##################################################

        if quick_test:
            # self.mom8_creator(image=img1+'.image.pbcor.tt0.subim', range_list=[0, 0.90])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")
        else:
            # self.mom8_creator(image=img1+'.image.pbcor.tt0.subim', range_list=[0, 0.90])
            print("skipping image creation - casa 6.1.3 hangs with imviewer")

        ##########################################
        # %% Generate Images [test_j1927_ql] end @
        ##########################################
        # not part of the jupyter scripts

        # save results for future analysis
        np.save(tstobj.id()+'.sourceflux.npy', curr_stats)
        np.save(tstobj.id()+'.beamstats.npy', beamstats_curr)

        # (f) Runtimes not significantly different relative to previous runs
        # don't test this in jupyter notebooks - runtimes differ too much between machines
        successr, reportr = tstobj.check_runtime(starttime, stats613['runtime'])
        report += reportr
        success = success and successr and self.th.check_final(report)

        tstobj.assertTrue(success, msg=report)

##############################################
##############################################

## List to be run0
def suite():
    return [test_j1302, test_j1927]

if __name__ == '__main__':
    unittest.main()
