"""
TestCases for casaxmlutil.py.

Unit test methods in this class are based on the tests of sdcal/sdfit.
If a method decorated by xml_constraints_injector of casaxmlutil.py are going to be added,
then you should consider adding some tests in this module for new constraints.
"""
import os
import shutil
import sys
import unittest

from casatasks import casalog, sdcal, sdfit
from casatasks.private.casaxmlutil import xml_constraints_injector
from casatools import ctsys

casalogpath = casalog.logfile()
testdir = 'cxutest'

sdcal_datapath = ctsys.resolve('unittest/sdcal')
sdcal_infiles = {'ps': 'uid___A002_X6218fb_X264.ms.sel', 'otf': 'uid___A002_X6218fb_X264.ms.sel.otfraster'}
sdcal_outfile = 'sdcal.out'
sdcal_applyfile = 'apply.cal'

sdfit_datapath = ctsys.resolve('unittest/sdfit')
sdfit_infiles = {'timebin': 'sdfit_tave.ms'}


class CasaxmlutilTest(unittest.TestCase):
    """Test class of casaxmlutil.py."""

    @classmethod
    def setUpClass(cls):
        """Create temoprary directory and copy files for tests."""
        cls.curdir = os.getcwd()
        if os.path.exists(testdir):
            shutil.rmtree(testdir)
        os.mkdir(testdir)
        os.chdir(testdir)
        for v in sdcal_infiles.values():
            if os.path.exists(v):
                shutil.rmtree(v)
            shutil.copytree(os.path.join(sdcal_datapath, v), v)
        for v in sdfit_infiles.values():
            if os.path.exists(v):
                shutil.rmtree(v)
            shutil.copytree(os.path.join(sdfit_datapath, v), v)

    @classmethod
    def tearDownClass(cls):
        """Reset current path and clear temporary files."""
        os.chdir(cls.curdir)
        if os.path.exists(testdir):
            shutil.rmtree(testdir)

    def tearDown(self):
        """Reset casalog and clear temporary files."""
        casalog.setlogfile(casalogpath)
        if os.path.exists(sdcal_outfile):
            shutil.rmtree(sdcal_outfile)
        if os.path.exists(sdcal_applyfile):
            shutil.rmtree(sdcal_applyfile)

    def __load_logfile(self, logfile):
        """Load the casalog file."""
        with open(logfile, 'r') as fp:
            for line in fp:
                yield line.rstrip()

    def __check_log(self, logfile, msg):
        """Check whether the casalog file contains the msg string or not."""
        for line in self.__load_logfile(logfile):
            if msg in line:
                return True
        return False

    def __test_positive(self, method, args: dict, desired: dict):
        """Execute a task with args, and check whether desired parameters have been overridden or not."""
        test_name = sys._getframe().f_back.f_code.co_name
        logfile = test_name + '.log'
        casalog.setlogfile(logfile)
        method(**args)
        for k, v in desired.items():
            self.assertTrue(self.__check_log(logfile, f"overrode argument: {k} -> '{v}'"))

    def __test_negative(self, method, args: dict):
        """Execute a task with args, and check desired parameters have been not overridden."""
        test_name = sys._getframe().f_back.f_code.co_name
        logfile = test_name + '.log'
        casalog.setlogfile(logfile)
        method(**args)
        self.assertFalse(self.__check_log(logfile, "overrode argument:"))

    @xml_constraints_injector
    def dummy(self, *args, **kwargs):
        """Raise ValueError when XML file loads."""
        return args, kwargs

    def test_dummy(self):
        """The method should raise error."""
        with self.assertRaises(ValueError):
            self.dummy()

    def test_sdcal_ps(self):
        """Test sdcal(calmode=ps)."""
        self.__test_negative(sdcal,
                             args={'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                                   'calmode': 'ps'})

    def test_sdcal_otfraster(self):
        """Test sdcal(calmode=otfraster)."""
        self.__test_positive(sdcal,
                             args={'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                                   'calmode': 'otfraster'},
                             desired={'intent': 'OBSERVE_TARGET#ON_SOURCE'})

    def test_sdcal_otfraster_not_override(self):
        """Test sdcal(calmode=otfraster, intent='..')."""
        self.__test_negative(sdcal,
                             args={'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                                   'calmode': 'otfraster', 'intent': 'OBSERVE_TARGET#ON_SOURCE'})

    def test_sdcal_otf(self):
        """Test sdcal(calmode=otf)."""
        self.__test_positive(sdcal,
                             args={'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                                   'calmode': 'otf'},
                             desired={'intent': 'OBSERVE_TARGET#ON_SOURCE'})

    def test_sdcal_otf_not_override(self):
        """Test sdcal(calmode=otf, intent='..')."""
        self.__test_negative(sdcal,
                             args={'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                                   'calmode': 'otf', 'intent': 'OBSERVE_TARGET#ON_SOURCE'})

    def test_sdcal_apply(self):
        """Test sdcal(calmode=apply)."""
        sdcal(infile=sdcal_infiles['ps'], calmode='tsys', outfile=sdcal_applyfile)
        self.__test_positive(sdcal,
                             args={'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                                   'calmode': 'apply', 'applytable': sdcal_applyfile},
                             desired={'interp': ''})

    def test_sdcal_apply_not_override(self):
        """Test sdcal(calmode=apply, interp='..')."""
        sdcal(infile=sdcal_infiles['ps'], calmode='tsys', outfile=sdcal_applyfile)
        self.__test_negative(sdcal,
                             args={'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                                   'calmode': 'apply', 'applytable': sdcal_applyfile, 'interp': 'nearest'})

    def test_sdcal_ps_apply(self):
        """Test sdcal(calmode=ps,apply)."""
        self.__test_positive(sdcal,
                             args={'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                                   'calmode': 'ps,apply'},
                             desired={'applytable': '', 'interp': ''})

    def test_sdcal_ps_apply_not_override(self):
        """Test sdcal(calmode=ps,apply) with applytabe and interp are not null."""
        sdcal(infile=sdcal_infiles['ps'], calmode='tsys', outfile=sdcal_applyfile)
        self.__test_negative(sdcal,
                             args={'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                                   'calmode': 'ps,apply', 'applytable': sdcal_applyfile,
                                   'interp': 'nearest'})

    def test_sdcal_tsys_apply(self):
        """Test sdcal(calmode=tsys,apply)."""
        self.__test_positive(sdcal,
                             args={'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                                   'calmode': 'tsys,apply'},
                             desired={'applytable': '', 'interp': ''})

    def test_sdcal_tsys_apply_not_override(self):
        """Test sdcal(calmode=tsys,apply) with applytabe and interp are not null."""
        sdcal(infile=sdcal_infiles['ps'], calmode='tsys', outfile=sdcal_applyfile)
        self.__test_negative(sdcal,
                             args={'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                                   'calmode': 'tsys,apply', 'applytable': sdcal_applyfile,
                                   'interp': 'nearest'})

    def test_sdcal_ps_tsys_apply(self):
        """Test sdcal(calmode=ps,tsys,apply)."""
        self.__test_positive(sdcal,
                             args={'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                                   'calmode': 'ps,tsys,apply'},
                             desired={'applytable': '', 'interp': ''})

    def test_sdcal_ps_tsys_apply_not_override(self):
        """Test sdcal(calmode=ps,tsys,apply) with applytabe and interp are not null."""
        sdcal(infile=sdcal_infiles['ps'], calmode='tsys', outfile=sdcal_applyfile)
        self.__test_negative(sdcal,
                             args={'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                                   'calmode': 'ps,tsys,apply', 'applytable': sdcal_applyfile,
                                   'interp': 'nearest'})

    def test_sdcal_otfraster_apply(self):
        """Test sdcal(calmode=otfraster,apply)."""
        self.__test_positive(sdcal,
                             args={'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                                   'calmode': 'otfraster,apply'},
                             desired={'applytable': '', 'interp': '', 'intent': 'OBSERVE_TARGET#ON_SOURCE'})

    def test_sdcal_otfraster_apply_not_override(self):
        """Test sdcal(calmode=otfraster,apply) with applytabe/interp/intent are not null."""
        sdcal(infile=sdcal_infiles['ps'], calmode='tsys', outfile=sdcal_applyfile)
        self.__test_negative(sdcal,
                             args={'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                                   'calmode': 'otfraster,apply', 'applytable': sdcal_applyfile,
                                   'interp': 'nearest', 'intent': 'OBSERVE_TARGET#ON_SOURCE'})

    def test_sdcal_otfraster_tsys_apply(self):
        """Test sdcal(calmode=otfraster,tsys,apply)."""
        self.__test_positive(sdcal,
                             args={'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                                   'calmode': 'otfraster,tsys,apply'},
                             desired={'applytable': '', 'interp': ''})

    def test_sdcal_otfraster_tsys_apply_not_override(self):
        """Test sdcal(calmode=otfraster,tsys,apply) with applytabe and interp are not null."""
        sdcal(infile=sdcal_infiles['ps'], calmode='tsys', outfile=sdcal_applyfile)
        self.__test_negative(sdcal,
                             args={'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                                   'calmode': 'otfraster,tsys,apply', 'applytable': sdcal_applyfile,
                                   'interp': 'nearest'})

    def test_sdcal_otf_tsys_apply(self):
        """Test sdcal(calmode=otf,tsys,apply)."""
        self.__test_positive(sdcal,
                             args={'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                                   'calmode': 'otf,tsys,apply'},
                             desired={'applytable': '', 'interp': ''})

    def test_sdcal_otf_tsys_apply_not_override(self):
        """Test sdcal(calmode=otf,tsys,apply) with applytabe and interp are not null."""
        sdcal(infile=sdcal_infiles['ps'], calmode='tsys', outfile=sdcal_applyfile)
        self.__test_negative(sdcal,
                             args={'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                                   'calmode': 'otf,tsys,apply', 'applytable': sdcal_applyfile,
                                   'interp': 'nearest'})

    def test_sdfit_timebin(self):
        """Test sdfit(timebin='1s'). If timebin is specified a value, then timespan is overridden by ''."""
        self.__test_positive(sdfit,
                             args={'infile': sdfit_infiles['timebin'],
                                   'datacolumn': 'float_data', 'nfit': [1], 'pol': 'XX',
                                   'timebin': '1s'},
                             desired={'timespan': ''})

    def test_sdfit_timebin_is_empty(self):
        """Test sdfit(timebin='')."""
        self.__test_negative(sdfit,
                             args={'infile': sdfit_infiles['timebin'],
                                   'datacolumn': 'float_data', 'nfit': [1], 'pol': 'XX',
                                   'timebin': ''})


if __name__ == '__main__':
    unittest.main()
