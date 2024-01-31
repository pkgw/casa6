#########################################################################
# test_task_gclean.py
#
# Copyright (C) 2018
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
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.imaging.tclean.html
#
# Each of the following categories (classes) has a set of tests within it.
#
#  test_onefield                 # basic tests, deconvolution algorithms
#  test_iterbot                   # iteration control options for mfs and cube
#  test_multifield               # multiple fields of same type and with different shapes/deconvolvers/gridders
#  test_stokes                    # multiple stokes planes, imaging with flagged correlations..
#  test_cube                      # all things cube. Spectral frame setup, handling empty channels, etc
#  test_widefield                # facets, wprojection, imagemosaic, mosaicft, awproject
#  test_mask                      # input mask options : regridding, mask file, automasking, etc
#  test_modelvis                # saving models (column/otf), using starting models, predict-only (setjy)
#  test_ephemeris                # ephemeris tests for gridder standard and mosaic, mode mfs and cubesource
#
#  To run the tests with python3 or casa in the command line
#
#  ./casa -c ./test_tclean.py
#  ./python3 ./test_tclean.py
#  ./casa -c ./casa6/casatestutils/runtest.py -h                                   # to see all options
#  ./casa -c ./casa6/casatestutils/runtest.py ./test_tclean.py                     # run the local test script
#  ./casa -c ./casa6/casatestutils/runtest.py ./test_tclean.py

# To run from within a casa 6 session:
#
#  from casatestutils import runtest
#  runtest.run(['test_tclean'])                                                 # pull test script from git trunk
#  runtest.run(['/path-to-test/test_tclean.py'])                                # run a local test script
#  runtest.run(['test_tclean[test_onefield_clark,test_onefield_hogbom]'])       # pull multiple test cases from git trunk
#  runtest.run(['test_tclean.py[test_onefield_clark,test_onefield_hogbom]'])    # run multiple test cases from local test script
#  See documentation for runtest.py in README.md of casatestutils
# To see the full list of tests :   grep "\"\"\" \[" test_tclean.py
#
#  These tests need data stored in casatestdata/unittest/tclean
#  The datasets are symliked to the above directory. If using cp to copy them locally,
#  Use cp -RH
#
#  For a developer build, to get the datasets locally
#
#  --- Get the test data repo :  svn co https://svn.cv.nrao.edu/svn/casatestdata casatestdata
#  --- Use ~/.casa/config.py to point to the casatestdata
# ########################################################################
# SKIPPED TESTS
# More tests were added to skip (as of 2019,04,26)
#
# (as of 2019.02.05 - Seven tests total)
# The following tests are currently skipped as the supports of the particular
# modes are not available in parallel mode yet
# =>
#     test_multifield_both_cube_diffshape
#     test_multifield_cube_mfs
#     test_multifield_cube_mtmfs
#
# The following tests in pricipal should be working but curently broken
# for parallel  until fixes to test or code are properly made.
# =>  test_multifield_facets_mfs
#     test_multifield_facets_mtmfs
#
# Added to skip at least for 5.5
#     test_cube_chanchunks
#     test_cube_chanchunks_savemodel (possible race conditions)
#     test_modelvis_2 (possible race conditions)
#     test_modelvis_3 (possible race conditions)
#     test_modelvis_5 (possible race conditions)
#     test_modelvis_6 (possible race conditions)
#     test_modelvis_7 (possible race conditions)
#     test_modelvis_8 (possible race conditions)
#     test_modelvis_9 (possible race conditions)
#     test_modelvis_10 (possible race conditions)
#     test_modelvis_11 (possible race conditions)
#     test_startmodel_with_mask_mfs(possible race conditions)
#     test_startmodel_with_mask_mtmfs(possible race conditions)

# Ressurected from skipping after some fixes
#     test_mask_5
#     test_iterbot_cube_2
#     test_multifield_both_cube
##########################################################################
#
#  Datasets
#
#  refim_twochan.ms : 2 channels, one 1Jy point source with spectral index of -1.0
#  refim_twopoints_twochan.ms : Two point sources, 1Jy and 5Jy, both with spectral index -1.0. For multifield tests.
#  refim_point.ms : 1-2 GHz, 20 channels, 1 spw, one 1Jy point source with spectral index -1.0.
#  refim_point_withline.ms : refim_point with a 'line' added into 3 channels (just topo)
#  refim_mawproject.ms : Two pointing wideband mosaic with 1 point source in between the two pointings
#  refim_mawproject_offcenter.ms : Two pointing wideband mosaic with 1 point source at center of one pointing
#  refim_point_stokes.ms : RR=1.0, LL=0.8, RL and LR are zero. Stokes I=0.9, V=0.1, U,Q=0.0
#  refim_point_linRL.ms : I=1, Q=2, U=3, V=4  in circular pol basis.
#  venus_ephem_test.ms : 7-point mosaic of Venus (ephemeris), Band 6, 1 spw, averaged to 1 chan
#
# List of test classes
#
# [test_onefield, test_iterbot, test_multifield,test_stokes, test_modelvis, test_cube, test_mask, test_startmodel, test_widefield,
# test_pbcor, test_mosaic_mtmfs, test_mosaic_cube, test_ephemeris, test_hetarray_imaging, test_wproject, test_errors_failures]
#
##########################################################################

import os
import sys
import shutil
import unittest
import inspect
import numpy as np
import operator

from casatools import ctsys, quanta, measures, image, vpmanager, calibrater
from casatasks import casalog, delmod, imsubimage, tclean, uvsub, imhead, imsmooth, immath, widebandpbcor, impbcor, flagdata, makemask
from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
from casatasks.private.imagerhelpers.parallel_imager_helper import PyParallelImagerHelper
from casatasks.private.imagerhelpers.summary_minor import SummaryMinor
from casatasks.private.imagerhelpers._gclean import gclean
from casatasks import impbcor, split, concat

from casatestutils.imagerhelpers import TestHelpers

_ia = image( )
_vp = vpmanager( )
_cb = calibrater( )
_qa = quanta( )
_me = measures( )

refdatapath = ctsys.resolve('unittest/tclean/')

defaultlogpath = casalog.logfile()

## Base Test class with Utility functions
class testref_base(unittest.TestCase):
     def setUp(self):
          self.epsilon = 0.05
          self.msfile = ""
          self.img = "tst"
          self.cfcache = 'cfcach'
          # To use subdir in the output image names in some tests (CAS-10937)
          self.img_subdir = 'refimager_tst_subdir'
          self.parallel = False
          self.nnode = 0
          if ParallelTaskHelper.isMPIEnabled():
              self.parallel = True
              self.PH = PyParallelImagerHelper()
              self.nnode = len(self.PH.getNodeList())

          self.th = TestHelpers()
          self.check_final = self.th.check_final

     def tearDown(self):
          """ don't delete it all """
          #self.delData()
          if casalog.logfile() != defaultlogpath:
              casalog.setlogfile(defaultlogpath)

     # Separate functions here, for special-case tests that need their own MS.
     def prepData(self,msname=""):
          os.system('rm -rf ' + self.img_subdir)
          os.system('rm -rf ' + self.img+'*')
          if msname != "":
               self.msfile=msname
          if (os.path.exists(self.msfile)):
               os.system('rm -rf ' + self.msfile)
          shutil.copytree(os.path.join(refdatapath,self.msfile), self.msfile)

     def prepCfcache(self,cfcache=""):
         if (os.path.exists(self.cfcache)):
               os.system('rm -rf ' + self.cfcache)
         if cfcache!="":
               self.cfcache=cfcache
         if (os.path.exists(self.cfcache)):
               os.system('rm -rf ' + self.cfcache)
         shutil.copytree(os.path.join(refdatapath,self.cfcache), self.cfcache)

     def delData(self,msname=""):
          if msname != "":
               self.msfile=msname
          if (os.path.exists(self.cfcache)):
               os.system('rm -rf ' + self.cfcache)
          if (os.path.exists(self.msfile)):
               os.system('rm -rf ' + self.msfile)
          os.system('rm -rf ' + self.img_subdir)
          os.system('rm -rf ' + self.img+'*')

     def prepInputmask(self,maskname=""):
          if maskname!="":
              self.maskname=maskname
          if (os.path.exists(self.maskname)):
              os.system('rm -rf ' + self.maskname)
          shutil.copytree(os.path.join(refdatapath,self.maskname), self.maskname)

     def prepInputTextFile(self, textfile=""):
          if textfile!="":
              self.textfile=textfile
              shutil.copy(os.path.join(refdatapath,self.textfile), self.textfile)


     def do_clean(self, vis, imagename, **kwargs):
         """
         Run _gclean in a manner similar to InteractiveClean() and test usage modes.
         Manual break at 20 iterations to prevent infinite loops.
         """
         clean = gclean(vis=vis, imagename=imagename, **kwargs)

         stopdesc, stopcode, majordone, nmajor, niter, retdict = clean.__next__()

         ncyc = 0
         while stopcode==0 or ncyc==0:
             stopdesc, stopcode, majordone, nmajor, niter, retdict = clean.__next__()

             ncyc = ncyc+1
             if ncyc == 20:
                 break ### This is just in case the stopping criteria fail in the test....

         return stopdesc, stopcode, majordone, nmajor, niter, retdict


     def fill_mask(self, maskname, fill_type='zero', channel=-1, stokes=-1):
         """
         Given an input mask, fill the mask with ones, zeros or a pre-defined
         box of ones. By default fill the mask in all channels and Stokes.

         Allowable fill_types : zeros, ones, box
         """

         if channel == -1:
             chan_slice = slice(None)
         else:
             chan_slice = slice(channel, channel+1)

         if stokes == -1:
             stokes_slice = slice(None)
         else:
             stokes_slice = slice(stokes, stokes+1)

         _ia.open(maskname)
         pix = _ia.getchunk()
         print("pix shape ", pix.shape)
         if fill_type == 'zero':
             pix[:, :, stokes_slice, chan_slice] = pix * 0.0
         elif fill_type == 'ones':
             pix[:, :, stokes_slice, chan_slice] = np.ones_like(pix)
         elif fill_type == 'box':
             pix[40:60,40:60, stokes_slice, chan_slice] = 1.0

             pix[0:40,0:40, stokes_slice, chan_slice] = 0.0
             pix[60:100,60:100, stokes_slice, chan_slice] = 0.0
             pix[0:40,60:100, stokes_slice, chan_slice] = 0.0
             pix[60:100,0:40, stokes_slice, chan_slice] = 0.0

         _ia.putchunk(pix)
         _ia.close()


     def calc_mask_sum(self, maskname):
         """
         Given an input mask, calculate the masksum
         """

         _ia.open(maskname)
         pix = _ia.getchunk()
         masksum = np.sum(pix)
         _ia.close()

         return masksum



class test_ic(testref_base):
    """
    Test iteration control options in gclean
    """

    def __init__(self, testref_base):
        super().__init__(testref_base)
        self.gclean = gclean

    # Test niter stopping criteria for cubes where niterdone > niter
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_niter_cube(self):
        """ [test_ic] Test_niter_cube : niter stopping criteria for cubes """

        self.prepData('refim_point.ms')

        stopdesc, stopcode, majordone, nmajor, niter, retdict = self.do_clean(vis=self.msfile, imagename=self.img, imsize=100, cell='10.0arcsec',
                                                                              specmode='cube', interpolation='nearest', nchan=5, start='1.0GHz', width='0.2GHz',
                                                                              pblimit=-1e-05, deconvolver='hogbom', niter=100, cycleniter=-1, cyclefactor=1, nmajor=3,
                                                                              threshold='0.01Jy', usemask='user', mask='')
        total_iterations = 0
        for nchan in range(5):
            # 'iterations' contains the cumulative sum, so first diff to get
            # the iterations per minor cycle Then sum to get the total number
            # of iterations.
            total_iterations += np.sum(np.diff(retdict['chan'][nchan][0]['iterations']))

        print("total iterations ", total_iterations)

        self.delData()

        # This should be the same as the number of major cycles done
        self.assertTrue(len(retdict['major']['cyclethreshold']) == 3)
        self.assertTrue(total_iterations == 156)
        self.assertTrue(stopcode == 1)


    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_nmajor(self):
        """ [test_ic] Test_nmajor : nmajor stopping criteria """

        self.prepData('refim_twochan.ms')

        stopdesc, stopcode, majordone, nmajor, niter, retdict = self.do_clean(vis=self.msfile, imagename=self.img, imsize=100, cell='10.0arcsec',
                                                                              specmode='mfs', interpolation='nearest', pblimit=-1e-05, deconvolver='hogbom',
                                                                              niter=100, cycleniter=10, cyclefactor=1, nmajor=3,
                                                                              threshold='0.01Jy', usemask='user', mask='')
        # This should be the same as the number of major cycles done
        # nmajor == 3 so it will trigger 4 major cycles total, including the
        # first one to make the initial residual image
        self.assertTrue(len(retdict['major']['cyclethreshold']) == 4)
        self.assertTrue(stopcode == 9)

        self.delData()


    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_niter_nomajor_nocycle(self):
        """ [test_ic] Test_niter_nomajor_nocycle : niter stopping criteria, nmajor=-1, cycleniter=-1 """

        self.prepData('refim_twochan.ms')

        stopdesc, stopcode, majordone, nmajor, niter, retdict = self.do_clean(vis=self.msfile, imagename=self.img, imsize=100, cell='10.0arcsec',
                                                                              specmode='mfs', interpolation='nearest', pblimit=-1e-05, deconvolver='hogbom',
                                                                              niter=100, cycleniter=-1, cyclefactor=1, nmajor=-1,
                                                                              threshold='0.01Jy', usemask='user', mask='')
        # This should be the same as the number of major cycles done
        self.assertTrue(len(retdict['major']['cyclethreshold']) == 4)
        self.assertTrue(stopcode == 1)


    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_niter_nocycle(self):
        """ [test_ic] Test_niter_nomajor_nocycle: niter stopping criteria, nmajor=2, cycleniter=-1  """

        self.prepData('refim_twochan.ms')

        stopdesc, stopcode, majordone, nmajor, niter, retdict = self.do_clean(vis=self.msfile, imagename=self.img, imsize=100, cell='10.0arcsec',
                                                                              specmode='mfs', interpolation='nearest', pblimit=-1e-05, deconvolver='hogbom',
                                                                              niter=100, cycleniter=-1, cyclefactor=1, nmajor=2,
                                                                              threshold='0.01Jy', usemask='user', mask='')
        # This should be nmajor+1
        self.assertTrue(len(retdict['major']['cyclethreshold']) == 3)
        self.assertTrue(stopcode == 9)

        self.delData()


    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_threshold(self):
        """ [test_ic] test_ic_threshold : threshold stopping criteria """

        self.prepData('refim_twochan.ms')

        stopdesc, stopcode, majordone, nmajor, niter, retdict = self.do_clean(vis=self.msfile, imagename=self.img, imsize=100, cell='10.0arcsec',
                                                                              specmode='mfs', interpolation='nearest', pblimit=-1e-05, deconvolver='hogbom',
                                                                              niter=100, cycleniter=-1, cyclefactor=1, nmajor=3,
                                                                              threshold='0.3Jy', usemask='user', mask='')

        self.assertTrue(len(retdict['major']['cyclethreshold']) == 2)
        self.assertTrue(stopcode == 2)
        self.delData()


    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_cube_cycleniter(self):
        """ [test_ic] test_ic_cube_cycleniter : cycleniter stopping criteria for cubes """

        self.prepData('refim_point.ms')

        stopdesc, stopcode, majordone, nmajor, niter, retdict = self.do_clean(vis=self.msfile, imagename=self.img, imsize=100, cell='10.0arcsec',
                                                                              specmode='cube', interpolation='nearest', nchan=5, start='1.0GHz', width='0.2GHz',
                                                                              pblimit=-1e-05, deconvolver='hogbom', niter=50, cycleniter=10, cyclefactor=1, nmajor=-1,
                                                                              threshold='0.01Jy', usemask='user', mask='')
        total_iterations = 0
        for nchan in range(5):
            # 'iterations' contains the cumulative sum, so first diff to get
            # the iterations per minor cycle Then sum to get the total number
            # of iterations.
            total_iterations += np.sum(np.diff(retdict['chan'][nchan][0]['iterations']))

        self.delData()

        # This should be the same as the number of major cycles done
        self.assertTrue(len(retdict['major']['cyclethreshold']) == 2)
        self.assertTrue(total_iterations == 50)
        self.assertTrue(stopcode == 1)




    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_mtmfs(self):
        """ [test_ic] test_ic_mtmfs : Check mtmfs naming and iteration control """

        self.prepData('refim_twochan.ms')

        stopdesc, stopcode, majordone, nmajor, niter, retdict = self.do_clean(vis=self.msfile, imagename=self.img, imsize=100, cell='10.0arcsec',
                                                                              specmode='mfs', interpolation='nearest', pblimit=-1e-05, deconvolver='mtmfs',
                                                                              niter=50, cycleniter=-1, cyclefactor=1, nmajor=-1,
                                                                              threshold='0.0Jy', usemask='user', mask='')

        self.delData()

        total_iterations = 0
        # 'iterations' contains the cumulative sum, so first diff to get
        # the iterations per minor cycle Then sum to get the total number
        # of iterations.
        total_iterations += np.sum(np.diff(retdict['chan'][0][0]['iterations']))

        self.assertTrue(len(retdict['major']['cyclethreshold']) == 4)
        self.assertTrue(stopcode == 1)
        self.assertTrue(total_iterations == 50)



    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_mfs_staticmask(self):
        """ [test_ic] test_ic_mfs_staticmask : Image MFS with a mask that does not change """

        self.prepData('refim_twochan.ms')


        stopdesc, stopcode, majordone, nmajor, niter, retdict = self.do_clean(vis=self.msfile, imagename=self.img, imsize=100,
                                                                               cell='10.0arcsec', specmode='mfs', interpolation='nearest',
                                                                               pblimit=-1e-05, deconvolver='hogbom', niter=50, cycleniter=-1,
                                                                               cyclefactor=1, nmajor=-1, threshold='0.0Jy', usemask='user',
                                                                               mask='circle[[50pix,50pix],10pix]')

        self.delData()

        total_iterations = 0
        # 'iterations' contains the cumulative sum, so first diff to get
        # the iterations per minor cycle Then sum to get the total number
        # of iterations.
        total_iterations += np.sum(np.diff(retdict['chan'][0][0]['iterations']))

        self.assertTrue(len(retdict['major']['cyclethreshold']) == 3)
        self.assertTrue(stopcode == 1)
        self.assertTrue(total_iterations == 50)


    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_automask(self):
        """ [test_ic] test_ic_automask : Image MFS with automasking """

        self.prepData('refim_twochan.ms')

        stopdesc, stopcode, majordone, nmajor, niter, retdict = self.do_clean(vis=self.msfile, imagename=self.img,
                                                                              imsize=100, cell='10.0arcsec', specmode='mfs',
                                                                              pblimit=-1e-05, deconvolver='hogbom', niter=100,
                                                                              cycleniter=10, cyclefactor=1, nmajor=3, threshold='0.0Jy',
                                                                              usemask='auto-multithresh', mask='')

        self.delData()

        total_iterations = 0
        total_iterations += np.sum(np.diff(retdict['chan'][0][0]['iterations']))

        self.assertTrue(len(retdict['major']['cyclethreshold']) == 4)
        self.assertTrue(stopcode == 9)
        self.assertTrue(total_iterations == 30)




    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_automask_staticmask(self):
        """ [test_ic] test_ic_automask_staticmask : Image MFS with automasking followed by a static mask """

        self.prepData('refim_twochan.ms')

        stopdesc1, stopcode1, majordone1, nmajor1, niter1, retdict1 = self.do_clean(vis=self.msfile, imagename=self.img,
                                                                              imsize=100, cell='10.0arcsec', specmode='mfs',
                                                                              pblimit=-1e-05, deconvolver='hogbom', niter=100,
                                                                              cycleniter=10, cyclefactor=1, nmajor=2, threshold='0.0Jy',
                                                                              usemask='auto-multithresh', mask='')
        masksum1 = self.calc_mask_sum(self.img+'.mask')

        # Now continue with a static mask
        stopdesc2, stopcode2, majordone2, nmajor2, niter2, retdict2 = self.do_clean(vis=self.msfile, imagename=self.img,
                                                                              imsize=100, cell='10.0arcsec', specmode='mfs',
                                                                              pblimit=-1e-05, deconvolver='hogbom', niter=100,
                                                                              cycleniter=10, cyclefactor=1, nmajor=2, threshold='0.0Jy',
                                                                              usemask='user', mask='')

        masksum2 = self.calc_mask_sum(self.img+'.mask')

        self.delData()

        self.assertTrue(masksum1 == masksum2)
        self.assertTrue(len(retdict1['major']['cyclethreshold']) == 3)
        self.assertTrue(len(retdict2['major']['cyclethreshold']) == 3)

        self.assertTrue(stopcode1 == 9)
        self.assertTrue(stopcode2 == 9)


    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_staticmask_automask(self):
        """ [test_ic] test_ic_staticmask_automask : Image MFS with static mask followed by automasking """

        self.prepData('refim_twochan.ms')

        stopdesc1, stopcode1, majordone1, nmajor1, niter1, retdict1 = self.do_clean(vis=self.msfile, imagename=self.img,
                                                                                imsize=100, cell='10.0arcsec', specmode='mfs',
                                                                                pblimit=-1e-05, deconvolver='hogbom', niter=100,
                                                                                cycleniter=10, cyclefactor=1, nmajor=2, threshold='0.0Jy',
                                                                                usemask='user', mask='circle[[50pix, 50pix], 2pix]')
        masksum1 = self.calc_mask_sum(self.img+'.mask')

        stopdesc2, stopcode2, majordone2, nmajor2, niter2, retdict2 = self.do_clean(vis=self.msfile, imagename=self.img,
                                                                                    imsize=100, cell='10.0arcsec', specmode='mfs',
                                                                                    pblimit=-1e-05, deconvolver='hogbom', niter=100,
                                                                                    cycleniter=10, cyclefactor=1, nmajor=2, threshold='0.0Jy',
                                                                                    usemask='auto-multithresh', mask='')

        masksum2 = self.calc_mask_sum(self.img+'.mask')
        self.delData()

        self.assertTrue(masksum1 == 13)
        self.assertTrue(masksum2 == 301)

        self.assertTrue(len(retdict1['major']['cyclethreshold']) == 3)
        self.assertTrue(len(retdict2['major']['cyclethreshold']) == 3)

        self.assertTrue(stopcode1 == 9)
        self.assertTrue(stopcode2 == 9)



    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_zeromask(self):
        """ [test_ic] test_ic_zeromask : Image MFS with a zero mask """

        self.prepData('refim_twochan.ms')


        stopdesc1, stopcode1, majordone1, nmajor1, niter1, retdict1 = self.do_clean(vis=self.msfile, imagename=self.img,
                                                                                    imsize=100, cell='10.0arcsec', specmode='mfs',
                                                                                    pblimit=-1e-05, deconvolver='hogbom', niter=100,
                                                                                    cycleniter=10, cyclefactor=1, nmajor=2, threshold='0.0Jy',
                                                                                    usemask='user', mask='')

        masksum1 = self.calc_mask_sum(self.img+'.mask')
        # Fill mask with zeros and start imaging again
        self.fill_mask(self.img+'.mask', fill_type='zero')


        stopdesc2, stopcode2, majordone2, nmajor2, niter2, retdict2 = self.do_clean(vis=self.msfile, imagename=self.img,
                                                                                    imsize=100, cell='10.0arcsec', specmode='mfs',
                                                                                    pblimit=-1e-05, deconvolver='hogbom', niter=100,
                                                                                    cycleniter=10, cyclefactor=1, nmajor=2, threshold='0.0Jy',
                                                                                    usemask='user', mask='')

        masksum2 = self.calc_mask_sum(self.img+'.mask')

        self.delData()


        self.assertTrue(masksum1 == 10000)
        self.assertTrue(masksum2 == 0)

        self.assertTrue(len(retdict1['major']['cyclethreshold']) == 3)
        self.assertTrue(len(retdict2['major']['cyclethreshold']) == 2)

        self.assertTrue(stopcode1 == 9)
        self.assertTrue(stopcode2 == 7)


    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_interaction(self):
        """ [test_ic] test_ic_interaction : Simulate interactively drawing and zeroing masks, and changing nmajor, niter and threshold """


        self.prepData('refim_twochan.ms')

        clean = self.gclean(vis=self.msfile, imagename=self.img, imsize=100, cell='10.0arcsec', specmode='mfs',
                            pblimit=-1e-05, deconvolver='hogbom', niter=100, cycleniter=10, cyclefactor=1, nmajor=2,
                            threshold='0.0Jy', usemask='user', mask='')

        # Trigger 1 + 2 major cycles, then confirm that the stopcode is 9
        # Fist major cycle just makes the residual
        # Make initial resiudal
        stopdesc, stopcode, majordone, nmajor, niter, retdict = clean.__next__()
        # Maj cycle 1
        stopdesc, stopcode, majordone, nmajor, niter, retdict = clean.__next__()
        # Maj cycle 2
        stopdesc, stopcode, majordone, nmajor, niter, retdict = clean.__next__()

        self.assertTrue(stopcode == 9)

        # Add more iterations and nmajor
        clean.update({'nmajor':9, 'niter':100, 'threshold':'0.01Jy'})

        # Zero out mask - should stop with stopcode = 7
        self.fill_mask(self.img+'.mask', fill_type='zero')
        stopdesc, stopcode, majordone, nmajor, niter, retdict = clean.__next__()

        self.assertTrue(stopcode == 7)

        # Fill mask with box
        self.fill_mask(self.img+'.mask', fill_type='box')
        # Update nmajor, niter and cycleniter
        clean.update({'nmajor':9, 'niter':100, 'threshold':'0.01Jy', 'cycleniter':50})
        stopdesc, stopcode, majordone, nmajor, niter, retdict = clean.__next__()
        stopdesc, stopcode, majordone, nmajor, niter, retdict = clean.__next__()
        stopdesc, stopcode, majordone, nmajor, niter, retdict = clean.__next__()

        self.assertTrue(stopcode == 1)

        # Set nmajor back to zero and it should stop with stopcode = 9
        clean.update({'nmajor':0, 'niter':100, 'threshold':'0.01Jy', 'cycleniter':50})
        stopdesc, stopcode, majordone, nmajor, niter, retdict = clean.__next__()

        self.assertTrue(stopcode == 9)

        self.delData()

        total_iterations = 0
        total_iterations += np.sum(np.diff(retdict['chan'][0][0]['iterations']))

        self.assertTrue(total_iterations == 120)


    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_cube_partialmask(self):
        """ [test_ic] test_ic_zeromask : Image cube with a single channel masked """

        self.prepData('refim_point.ms')

        stopdesc1, stopcode1, majordone1, nmajor1, niter1, retdict1 = self.do_clean(vis=self.msfile, imagename=self.img, imsize=100, cell='10.0arcsec', specmode='cube',
                                                                              interpolation='nearest', nchan=5, start='1.0GHz', width='0.2GHz', pblimit=-1e-05,
                                                                              deconvolver='hogbom', niter=10, cycleniter=-1, cyclefactor=1, nmajor=3,
                                                                              threshold='0.01Jy', usemask='user', mask='')

        masksum1 = self.calc_mask_sum(self.img+'.mask')
        # Fill mask with zeros and start imaging again
        self.fill_mask(self.img+'.mask', fill_type='box', channel=1)

        # Clean down until single channel reaches threshold
        stopdesc2, stopcode2, majordone2, nmajor2, niter2, retdict2 = self.do_clean(vis=self.msfile, imagename=self.img,
                                                                                    imsize=100, cell='10.0arcsec', specmode='cube',
                                                                                    interpolation = 'nearest', nchan=5, start='1.0GHz', width='0.2GHz',
                                                                                    pblimit=-1e-05, deconvolver='hogbom', niter=100,
                                                                                    cycleniter=10, cyclefactor=1, nmajor=2, threshold='1.0Jy',
                                                                                    usemask='user', mask='')

        masksum2 = self.calc_mask_sum(self.img+'.mask')
        self.fill_mask(self.img+'.mask', fill_type='box', channel=2)
        stopdesc3, stopcode3, majordone3, nmajor3, niter3, retdict3 = self.do_clean(vis=self.msfile, imagename=self.img,
                                                                                    imsize=100, cell='10.0arcsec', specmode='cube',
                                                                                    interpolation = 'nearest', nchan=5, start='1.0GHz', width='0.2GHz',
                                                                                    pblimit=-1e-05, deconvolver='hogbom', niter=100,
                                                                                    cycleniter=10, cyclefactor=1, nmajor=2, threshold='1.0Jy',
                                                                                    usemask='user', mask='')
        masksum3 = self.calc_mask_sum(self.img+'.mask')

        self.delData()

        self.assertTrue(masksum1 == 0)
        self.assertTrue(masksum2 == 400)
        self.assertTrue(masksum3 == 800)

        self.assertTrue(len(retdict1['major']['cyclethreshold']) == 2)
        self.assertTrue(len(retdict2['major']['cyclethreshold']) == 2)
        self.assertTrue(len(retdict3['major']['cyclethreshold']) == 2)

        self.assertTrue(stopcode1 == 7)
        self.assertTrue(stopcode2 == 2)
        self.assertTrue(stopcode3 == 2)


    # Test niter stopping criteria for cubes where niterdone > niter
    @unittest.skipIf(ParallelTaskHelper.isMPIEnabled(), "gclean doesn't work with mpi")
    def test_ic_niter_cube(self):
        """ [test_ic] Test_niter_cube : niter stopping criteria for cubes """

        total_iterations = 0
        for nchan in range(5):
            # 'iterations' contains the cumulative sum, so first diff to get
            # the iterations per minor cycle Then sum to get the total number
            # of iterations.
            total_iterations += np.sum(np.diff(retdict['chan'][nchan][0]['iterations']))

        print("total iterations ", total_iterations)

        self.delData()

        # This should be the same as the number of major cycles done
        self.assertTrue(len(retdict['major']['cyclethreshold']) == 3)
        self.assertTrue(total_iterations == 156)
        self.assertTrue(stopcode == 1)


