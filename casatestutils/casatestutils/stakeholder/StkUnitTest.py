import os
import copy
import shutil
from datetime import datetime
from collections import Iterable
import numpy as np
import unittest
import glob

from casatools import ctsys, table
from casatools import image as _ia
from casatasks import casalog, tclean
from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper
from casatestutils.imagerhelpers import TestHelpers

class StkUnitTest(unittest.TestCase):
    """ Adds some stakeholder test specific extensions to the general unit test class """
    
    def setUp(self):
        super().setUp()
        self.imgs = []
        self.th = TestHelpers()
        self.tb = table()
        self.ia = _ia()
        self.vis = "" # measurement set name
        self.parallel = ParallelTaskHelper.isMPIEnabled()
        self._clean_imgs_exist_dict()

    def tearDown(self):
        super().tearDown()
        if os.getenv("CACHE_PARTIAL_RESULTS") != "true":
            self.delData()

    def delData(self):
        """ Clean up generated data. """
        del_files = []
        if self.vis != "" and os.path.exists(self.vis):
            del_files.append(self.vis)
        for img in self.imgs:
            img_files = glob.glob(img+'*')
            del_files += img_files
        for f in del_files:
            shutil.rmtree(f)

    def _copy_file_or_dir(self, src, dst):
        if (os.path.isdir(src)):
            shutil.copytree(src, dst)
        else:
            shutil.copy2(src, dst)

    def prepData(self, msname, data_path_dir, *copyargs):
        """ Copies the given measurement set (and other copyargs) to the current directory.

        This function fulfills a different purpose in the case that the environment variable
          USE_PARTIAL_RESULTS == "true"
        In that case, instead of copying the ms and copyargs, it does the equivalent of
          cp -rp ../../partial_results/* ./
        This can be useful to, for example, run tclean once, then quickly iterate on the
        analysis portion of the tests.

        Args:
          msname: name of the primary measurement set
          data_path_dir: subdirectory containing the measurement set and any other copyargs
          copyargs: additional directories to be copied from data_path_dir to the current directory
          USE_PARTIAL_RESULTS: env var, "true" to copy partial results
        """
        self.vis = msname

        if os.getenv("USE_PARTIAL_RESULTS") != "true":
            # clean run
            data_path_dir = ctsys.resolve(data_path_dir)
            mssrc = os.path.join(data_path_dir, self.vis)
            casalog.post(f"{mssrc} => {self.vis}", "INFO")
            shutil.copytree(mssrc, self.vis)

            for copydir in copyargs:
                copysrc = os.path.join(data_path_dir, copydir)
                casalog.post(f"{copysrc} => {copydir}", "INFO")
                shutil.copytree(copysrc, copydir)
        else:
            # continue running with partially computed results (eg, ran tclean last time, now check the values)
            from os.path import dirname, join
            import sys
            fromdir = join( dirname(dirname(os.getcwd())), "partial_results" )
            skipfiles = ["__pycache__"]
            files = list(filter(lambda x: x not in skipfiles, os.listdir(fromdir)))
            casalog.post(f"Restorting partial results [{len(files)}]", "INFO")
            for i in range(len(files)):
                casalog.post(f"{i}: {files[i]}", "INFO")
                self._copy_file_or_dir(join(fromdir, files[i]), files[i])

    def check_img_exists(self, img):
        """ Returns true if the image exists. A report is collected internally, to be returned as a group report in get_imgs_exist_results(...).

        Args:
          img: The name of the image directory to check for

        Returns:
          success bool

        See also:
          _clean_imgs_exist_dict()
          get_imgs_exist_results()
        """
        exists = self.th.image_exists(img)
        success, report = self.th.check_val(exists, True, valname=f"image_exists('{img}')", exact=True, testname=self._testMethodName)
        if not exists:
            # log immediately: missing images could cause the rest of the test to fail
            casalog.post(report, "SEVERE")
        self.imgs_exist['successes'].append(success)
        self.imgs_exist['reports'].append(report)
        return success

    def get_imgs_exist_results(self):
        """ Get a single collective result of check_img_exists(...)

        Returns:
          (success bool, report string)
        """
        success = all(self.imgs_exist['successes'])
        report = "".join(self.imgs_exist['reports'])
        return success, report

    def _clean_imgs_exist_dict(self):
        """ Clean the arrays that hold onto the list of images to check for existance.

        See also: check_img_exists(...)
        """
        self.imgs_exist = { 'successes':[], 'reports':[] }

    def _nparray_to_list(self, val):
        if isinstance(val, np.ndarray) and val.ndim > 0:
            val = list(val)
            for i in range(len(val)):
                val[i] = self._nparray_to_list(val[i])
        return val

    def check_diff(self, actual, expected, diff, valname, desired_diff, max_diff, rms=None):
        """ Check that the given difference is within tolerance.

        Compare the given difference to the desired and maximum tolerances.
        Logs a warning if outside of desired bounds.

        Args:
          actual: measured value(s), to be included in the report
          expected: truth value(s), to be included in the report
          diff: the difference between the actual and expected value(s)
          valname: name of the value being compared, to be included in the report
          desired_diff: the tolerance for expected we'd like to stay within, but that is not strictly required
          max_diff: the required tolerance for expected that causes this test to fail
          rms: the root mean square of the noise floor, to be included in the report

        Returns:
          (success bool, report string)
        """

        # only worry about comparing the maximum value
        val = diff
        if isinstance(diff, Iterable):
            val = max(diff)
        
        # convert numpy arrays to lists so that the logs get printed on a single line
        actual = self._nparray_to_list(actual)
        expected = self._nparray_to_list(expected)
        diff = self._nparray_to_list(diff)

        # get some values
        success = np.all(val <= max_diff)
        testname = self._testMethodName
        correctval = f"< {max_diff}"

        # generate the report
        rms_str = "" if (rms == None) else f", rms: {rms}"
        if (desired_diff != None) and (val > desired_diff):
            casalog.post(f"Warning, {valname}: {diff} vs desired {desired_diff}, (actual: {actual}, expected: {expected}{rms_str})", "WARN")
        report = "[ {} ] {} is {} ( {} : should be {}{})\n".format(testname, valname, str(diff), self.th.verdict(success), str(correctval), rms_str )
        report = report.rstrip() + f" (raw actual/expected values: {actual}/{expected})\n"
        if not success:
            casalog.post(report, "WARN") # easier to read this way than in an assert statement

        return success, report

    def check_fracdiff(self, actual, expected, valname, desired_diff=0.05, max_diff=0.1):
        """ Calls check_diff( diff=(actual-expected)/expected )

        Returns:
          (success bool, report string)
        """
        fracdiff=abs(actual-expected)/abs(expected)
        return self.check_diff(actual, expected, fracdiff, valname, desired_diff, max_diff)

    def check_metrics_flux(self, measured, expected, valname, rms_or_std, desired_diff=0.05, max_diff=0.1, nsigma=2):
        """ Logs a warning if outside of desired bounds, returns False if outside required bounds
        
        Check that the given value(s) are within a reasonable tolerance of the truth value(s),
        while taking the rms noise floor into account. The default tolerances are:
        5% desired, 10% required, as from https://drive.google.com/file/d/1zw6UeDEoXoxM05oFg3rir0hrCMEJMxkH/view and https://open-confluence.nrao.edu/display/VLASS/Updated+VLASS+survey+science+requirements+and+parameters

        Args:
          a_meas: measured value(s)
          a_true: truth value(s)
          valname: name of the value to be printed in the report
          rms_or_std: root mean square of the noise floor
          desired_diff: the tolerance for a_true we'd like to stay within, but that is not strictly required
          max_diff: the required tolerance for a_true that causes this test to fail
          nsigma: within how many standard deviations of the noise floor is a value considered noise

        Returns:
          (success bool, report string)

        See also:
          https://casadocs.readthedocs.io/en/latest/notebooks/synthesis_imaging.html#Options-in-CASA-for-wideband-imaging --> Search for "Calculating Error in Spectral Index"
          Eqn 39 of https://www.aanda.org/index.php?option=com_article&access=doi&doi=10.1051/0004-6361/201117104&Itemid=129#S29
        """
        std = rms_or_std

        #################################################
        ###  [A] Tolerance on the metric of relative error : Based only on image noise levels
        #################################################
        calc = np.sign(expected) * ( np.abs(expected) + nsigma*std )
        rel_error = np.abs(( calc - expected ))/ np.maximum(np.abs(expected),nsigma*std)

        #################################################
        ### [B] Empirical tolerances
        #################################################
        emp_error = max_diff ## 0.1 = 10% relative error

        #################################################
        ### Tolerance to use for the tests
        #################################################
        tol_flux = np.maximum(rel_error, emp_error)  

        #################################################
        ### Calculate the metrics for the input measured values
        #################################################
        rel_error_flux = np.abs( ( measured - expected) ) / np.maximum(np.abs(expected),nsigma*std)
        
        return self.check_diff(measured, expected, rel_error_flux, valname, desired_diff, tol_flux)

    def check_metrics_alpha(self, a_meas, a_true, valname, rmss_or_stds, desired_diff=0.1, max_diff=0.2, nsigma=2):
        """ Logs a warning if outside of desired bounds, returns False if outside required bounds
        
        Check that the given value(s) are within a reasonable tolerance of the truth value(s),
        while taking the rms noise floor into account. The default tolerances are:
        0.1 desired, 0.2 required, as from https://drive.google.com/file/d/1zw6UeDEoXoxM05oFg3rir0hrCMEJMxkH/view and https://open-confluence.nrao.edu/display/VLASS/Updated+VLASS+survey+science+requirements+and+parameters

        Args:
          a_meas: measured value(s)
          a_true: truth value(s)
          valname: name of the value to be printed in the report
          rmss_or_stds: rms or [rms, rms], root mean square of the noise floor for the observed .image or .image.tt0, .image.tt1
          desired_diff: the tolerance for a_true we'd like to stay within, but that is not strictly required
          max_diff: the required tolerance for a_true that causes this test to fail
          nsigma: within how many standard deviations of the noise floor is a value considered noise

        Returns:
          (success bool, report string)

        See also:
          https://casadocs.readthedocs.io/en/latest/notebooks/synthesis_imaging.html#Options-in-CASA-for-wideband-imaging --> Search for "Calculating Error in Spectral Index"
          Eqn 39 of https://www.aanda.org/index.php?option=com_article&access=doi&doi=10.1051/0004-6361/201117104&Itemid=129#S29
        """
        if type(rmss_or_stds) != list:
            rmss_or_stds = [rmss_or_stds, rmss_or_stds]
        std0, std1 = rmss_or_stds[0], rmss_or_stds[1]

        #################################################
        ###  [A] Tolerance on the metric of relative error : Based only on image noise levels
        #################################################
        #       dIa =             Ia   *    sqrt( (    dI0  /  I0    )**2 + (    dI1  /  I1    )**2 )
        a_error_std = np.abs(   a_true * np.sqrt( (nsigma*std0/a_true)**2 + (nsigma*std1/a_true)**2 )   )

        #################################################
        ### [B] Empirical tolerances
        #################################################
        emp_a = max_diff ## probably 0.2, absolute error in spectral index 

        #################################################
        ### Tolerance to use for the tests
        #################################################
        tol_a =np.max([a_error_std, emp_a])

        #################################################
        ### Calculate the metrics for the input measured values
        #################################################
        abs_error_a = np.abs( a_meas - a_true )
        
        return self.check_diff(a_meas, a_true, abs_error_a, valname, desired_diff, tol_a)

    def check_metrics_alpha_fitted(self, a_meas, a_true, valname, pcov, desired_diff=0.1, max_diff=0.2, nsigma=2):
        """ For checking alpha, accounting for errors in fitting across several spws with scipy.optimize.curve_fit

        For example:
          from scipy.optimize import curve_fit

          freqs = [ 2.028, 2.796, 3.594 ]
          fluxes = [
              [ imstat(imagename=f"spw{i}.image.tt0",box=box,stokes=s) for s in "IQUV" ],
              for i in range(3)
            ]
          logfreqs=np.log10(freqs)
          logfluxes=np.log10(fluxes)
  
          def func(x, a, b):
             nu_0=3.0
             return a*(x-np.log10(nu_0))+b
          popt, pcov = curve_fit(func, logfreqs, logfluxes)

          alpha = popt[0]
          f_nu0 = 10**popt[1]
          success, report = tstobj.check_metrics_alpha_fitted(a_meas=alpha, a_true=0.3, valname="alpha", pcov=pcov)

        Returns:
          (success bool, report string)
        """
        std_a = nsigma * np.sqrt(np.diag(pcov)[0])
        emp_a = max_diff ## probably 0.2, absolute error in spectral index 
        tol_a = np.max([std_a, emp_a])
        abs_error_a = np.abs( a_meas - a_true )

        return self.check_diff(a_meas, a_true, abs_error_a, valname, desired_diff, tol_a)

    def check_column_exists(self, colname):
        """ Verifies that the given column exists in the self.vis measurement set.

        Returns:
          (success bool, report string)"""
        self.tb.open(self.vis)
        cnt = self.tb.colnames().count(colname)
        self.tb.done()
        self.tb.close()
        success, report = self.th.check_val(cnt, 1, valname=f"count('{colname}')", exact=True, testname=self._testMethodName)
        if not success:
            casalog.post(report, "WARN") # easier to read this way than in an assert statement
        return success, report

    def check_runtime(self, starttime, exp_runtime):
        """ Verifies that the runtime is within 10% of the expected runtime.

        Probably only valid when running on the same hardware as was used to measure the previous runtime.

        Returns:
          (success bool, report string)
        """
        endtime         = datetime.now()
        runtime         = (endtime-starttime).total_seconds()
        success, report = self.th.check_val(runtime, exp_runtime, valname="runtime", exact=False, epsilon=0.1, testname=self._testMethodName)
        if not success:
            casalog.post(report, "WARN") # easier to read this way than in an assert statement
        return success, report

    def get_params_as_dict(self, **wargs):
        """ Get the parameters called of a given function as a dictionary of parameter name to value.

        This can be useful to, for example, get a dictionary of parameter values
        for a previous call to tclean from the casalogs:
          self.get_params_as_dict(vis='J1927_12fields.ms', selectdata=True, field='', spw='', timerange='', uvrange='', antenna='', scan='', observation='', intent='OBSERVE_TARGET#UNSPECIFIED', datacolumn='data', imagename='VLASS1.2.ql.T26t15.J1927.10.2048.v1.I.iter0', imsize=[7290, 7290], cell='1.0arcsec', phasecenter='19:27:30.443 +61.17.32.898', stokes='I', projection='SIN', startmodel='', specmode='mfs', reffreq='3.0GHz', nchan=-1, start='', width='', outframe='LSRK', veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=False, gridder='mosaic', facets=1, psfphasecenter='', chanchunks=1, wprojplanes=1, vptable='', mosweight=False, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='', usepointing=False, computepastep=360.0, rotatepastep=360.0, pointingoffsetsigdev=[], pblimit=0.2, normtype='flatnoise', deconvolver='mtmfs', scales=[0], nterms=2, smallscalebias=0.0, restoration=False, restoringbeam='common', pbcor=False, outlierfile='', weighting='briggs', robust=1.0, noise='1.0Jy', npixels=0, uvtaper=[], niter=0, gain=0.1, threshold='0.0mJy', nsigma=0.0, cycleniter=-1, cyclefactor=1.0, minpsffraction=0.05, maxpsffraction=0.8, interactive=0, usemask='user', mask='', pbmask=0.0, sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0, smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True, minpercentchange=-1.0, verbose=False, fastnoise=True, restart=True, savemodel='none', calcres=True, calcpsf=True, parallel=False)
        """
        return dict(wargs)

    def print_task_diff_params(self, fname, act_pars : dict, exp_pars : dict):
        """ Compare the parameter values for the "act_pars" actual parameters
        and the "exp_pars" expected parameters. Print the parameters
        that are different and what their actual/expected values are.

        This can be useful to, for example, verify that all the parameter values
        are the same between a tclean call from a stakeholder's script file and
        a tclean call from a test script.
        """
        same_par_vals = []
        diff_par_vals = []
        diff_par_strs = []
        new_par_vals = []
        new_par_strs = []

        for pname in act_pars:
            par_found = False
            aval_differs = True
            aval = act_pars[pname]
            aval_str = f"'{aval}'" if (type(aval) == str) else str(aval)
            xval = None if pname not in exp_pars else exp_pars[pname]
            xval_str = f"'{xval}'" if (type(xval) == str) else str(xval)
            if pname in exp_pars:
                par_found = True
                if aval == exp_pars[pname]:
                    same_par_vals.append(pname)
                    aval_differs = False
            if not par_found:
                new_par_vals.append(pname)
                new_par_strs.append(f"{pname}={aval_str}")
            elif aval_differs:
                diff_par_vals.append(pname)
                diff_par_strs.append(f"{pname}={aval_str}/{xval_str}")

        diff_pars_str = ", ".join(diff_par_strs)
        new_pars_str = ", ".join(new_par_strs)

        casalog.post(    f"These parameters are different/new: {diff_par_vals+new_par_vals}", "INFO")
        if len(diff_pars_str) > 0:
            casalog.post(f"                 (actual/expected): {diff_pars_str}", "INFO")
        if len(new_par_vals) > 0:
            casalog.post(f"                          new pars: {new_pars_str}", "INFO")

    def _run_tclean(self, **wargs):
        """ Tracks the "imagename" in self.imgs (for cleanup), checks for mask existance, and runs tclean.

        If the env var USE_PARTIAL_RESULTS == "true", then don't run tclean.
        """
        if ('imagename' in wargs):
            img = wargs['imagename']
            if (img not in self.imgs):
                self.imgs.append(img)
        if ('mask' in wargs) and (wargs['mask'] != ''):
            if not os.path.exists(wargs['mask']):
                raise RuntimeError(f"Error: trying to run tclean with nonexistant mask {wargs['mask']}")
        try:
            if os.getenv("USE_PARTIAL_RESULTS") != "true":
                return tclean(**wargs)
                pass
        except:
            # self.print_tclean(**wargs)
            raise

    def run_tclean(self, vis='', selectdata=True, field='', spw='', timerange='', uvrange='', antenna='',
                    scan='', observation='', intent='', datacolumn='corrected', imagename='', imsize=[100],
                    cell=['1arcsec'], phasecenter='', stokes='I', projection='SIN', startmodel='',
                    specmode='mfs', reffreq='', nchan=- 1, start='', width='', outframe='LSRK',
                    veltype='radio', restfreq=[], interpolation='linear', perchanweightdensity=True,
                    gridder='standard', facets=1, psfphasecenter='', wprojplanes=1, vptable='',
                    mosweight=True, aterm=True, psterm=False, wbawp=True, conjbeams=False, cfcache='',
                    usepointing=False, computepastep=360.0, rotatepastep=360.0, pointingoffsetsigdev=[],
                    pblimit=0.2, normtype='flatnoise', deconvolver='hogbom', scales='', nterms=2,
                    smallscalebias=0.0, restoration=True, restoringbeam='', pbcor=False, outlierfile='',
                    weighting='natural', robust=0.5, noise='1.0Jy', npixels=0, uvtaper=[], niter=0,
                    gain=0.1, threshold=0.0, nsigma=0.0, cycleniter=- 1, cyclefactor=1.0, minpsffraction=0.05,
                    maxpsffraction=0.8, interactive=False, usemask='user', mask='', pbmask=0.0,
                    sidelobethreshold=3.0, noisethreshold=5.0, lownoisethreshold=1.5, negativethreshold=0.0,
                    smoothfactor=1.0, minbeamfrac=0.3, cutthreshold=0.01, growiterations=75, dogrowprune=True,
                    minpercentchange=- 1.0, verbose=False, fastnoise=True, restart=True, savemodel='none',
                    calcres=True, calcpsf=True, psfcutoff=0.35, parallel=None, compare_tclean_pars=None):
        """ Runs tclean with the default parameters from v6.4.0
        If the 'compare_tclean_pars' dict is provided, then compare these values to the other parameters of this function.

        See also:
          self.run_tclean(...)
        """
        parallel = (self.parallel) if (parallel == None) else (parallel)
        run_tclean_pars = locals()
        run_tclean_pars = {k:run_tclean_pars[k] for k in filter(lambda x: x not in ['self', 'compare_tclean_pars', 'psfcutoff'] and '__' not in x, run_tclean_pars.keys())}
        if (compare_tclean_pars != None):
            self.print_task_diff_params('run_tclean', act_pars=run_tclean_pars, exp_pars=compare_tclean_pars)
        return self._run_tclean(**run_tclean_pars)

    def resize_mask(self, maskname, outputmaskname, shape=[2000,2000,1,1]):
        """ Resizes a .image mask from its current shape to the given shape """
        from skimage.transform import resize # from package scikit-image

        # scrub the input
        # a trailing '/' on the image name causes calc to give weird results
        maskname = maskname.rstrip(" \t\n\r/\\")
        if (len(shape) < 4): # ra, dec, chan, pol
            shape += [1 for i in range(4-len(shape))]
        if shape[2] != 1 or shape[3] != 1:
            raise RuntimeError("Error: image must have length 1 in the third (chan) and fourth (pol) dimensions")

        # get the shape
        self.ia.open(maskname)
        try:
            inshape = self.ia.shape()
            pixeltype = self.ia.pixeltype()
            inpixels = self.ia.getregion()
        finally:
            self.ia.close()
            self.ia.done()

        # populate some pixels
        pixels = resize(inpixels, shape)
        for r in range(shape[0]):
            for d in range(shape[1]):
                pixels[r][d][0][0] = 0 if (pixels[r][d][0][0] < 0.5) else 1

        # create the new outputmask
        if (pixeltype == 'dcomplex'):
            pixeltype = 'cd'
        else:
            pixeltype = pixeltype[0] # for 'f'loat, 'd'ouble, or 'c'omplex
        if os.path.exists(outputmaskname):
            shutil.rmtree(outputmaskname)
        self.ia.fromarray(outputmaskname, pixels=pixels, type=pixeltype)
        self.ia.close()
        self.ia.done()