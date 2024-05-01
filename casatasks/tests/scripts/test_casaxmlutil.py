"""
TestCases for casaxmlutil.py.

Unit test methods in this class are based on the tests of sdcal/sdfit.
If you add a new task decorated by xml_constraints_injector or decorate an existing task
with xml_constraints_injector of casaxmlutil.py, then you should consider adding some tests
in this module for new constraints.
"""
import inspect
import os
import shutil
from typing import NamedTuple, Callable
import unittest

from casatasks import casalog, sdcal, sdfit
from casatasks.private.casaxmlutil import xml_constraints_injector
from casatools import ctsys

casalogpath = casalog.logfile()
testdir = 'cxutest'
SUCCESS = True
FAIL = False


class Testdata(NamedTuple):
    """Testdatapath class."""

    datapath: str
    infiles: dict
    tempfiles: dict


class CasaxmlutilTest(unittest.TestCase):
    """Test class of casaxmlutil.py."""

    _DUMMY = 'dummy'
    # testdata must be defined
    testdata = Testdata(_DUMMY, {}, {})
    task = None

    @classmethod
    def setUpClass(cls):
        """Create temporary directory and copy files for tests."""
        assert cls.testdata.datapath != cls._DUMMY, 'datapath was not properly configured'
        assert cls.task is not None, 'task was not specified'

        cls.curdir = os.getcwd()
        if os.path.exists(testdir):
            shutil.rmtree(testdir)
        os.mkdir(testdir)
        os.chdir(testdir)
        for td_path in cls.testdata.infiles.values():
            if os.path.exists(td_path):
                shutil.rmtree(td_path)
            shutil.copytree(os.path.join(cls.testdata.datapath, td_path), td_path)

    @classmethod
    def tearDownClass(cls):
        """Reset current path and clear temporary files."""
        os.chdir(cls.curdir)
        if os.path.exists(testdir):
            shutil.rmtree(testdir)

    def tearDown(self):
        """Reset casalog and clear temporary files."""
        casalog.setlogfile(casalogpath)
        for tempfile in self.testdata.tempfiles.values():
            if os.path.exists(tempfile):
                shutil.rmtree(tempfile)

    def check_log(self, logfile: str, msg: str):
        """Check whether the casalog file contains the msg string or not."""
        with open(logfile, 'r') as fp:
            for _ in filter(lambda x: msg in x, fp):
                return True
        return False

    def __test(self, method: Callable, args: dict):
        """Execute a task with args and return logfile name."""
        logfile = inspect.stack()[3].function + '.log'
        casalog.setlogfile(logfile)
        method(**args)
        self.assertTrue(os.access(logfile, os.R_OK))
        return logfile

    def positive_test(self, method: Callable, args: dict, desired: dict):
        """Execute a task with args, and check whether desired parameters have been overridden or not."""
        logfile = self.__test(method, args)
        for k, v in desired.items():
            self.assertTrue(self.check_log(logfile, f"overrode argument: {k} -> '{v}'"))

    def negative_test(self, method: Callable, args: dict):
        """Execute a task with args, and check desired parameters have been not overridden."""
        logfile = self.__test(method, args)
        self.assertFalse(self.check_log(logfile, "overrode argument:"))


class DummyGoImpClass:
    """DummyClass for FundamentalTest."""

    def __call__(self, *args, **kwargs):
        """Execute as an instance itself."""
        return self.dummy()

    @xml_constraints_injector
    def dummy(self, *args: list, **kwargs: dict):
        """Return True only."""
        return True


class RaiseExceptionClass:
    """TestClass for Exception."""

    def __call__(self, *args, **kwargs):
        """Execute as an instance itself."""
        return self.dummy()

    @xml_constraints_injector
    def dummy(self, *args: list, **kwargs: dict):
        """Return True only."""
        raise CasaXMLException


class CasaXMLException(BaseException):
    """Dummy Exception class."""


class FundamentalTest(unittest.TestCase):
    """Test fundamental use of xml_constraints_injector."""

    def test_dummy_pass(self):
        """The method should pass."""
        _tmp = DummyGoImpClass()
        self.assertTrue(_tmp())

    def test_raise_exception(self):
        """Check propagation of Exception"""
        _tmp = RaiseExceptionClass()
        with self.assertRaises(CasaXMLException):
            _tmp()


class SDCalTest(CasaxmlutilTest):
    """Test the constraints of sdcal."""

    testdata = Testdata(ctsys.resolve('unittest/sdcal'),
                        {'ps': 'uid___A002_X6218fb_X264.ms.sel',
                         'otf': 'uid___A002_X6218fb_X264.ms.sel.otfraster'},
                        {'outfile': 'sdcal.out',
                         'applyfile': 'apply.cal'})
    task = sdcal
    success = SUCCESS
    fail = FAIL

    def __test_sdcal(self, whether: bool, args: dict, desired: dict=None, prepare_data: bool=False):
        """Test sdcal with parameters, prepare input data if needed."""
        if prepare_data:
            self.task(infile=self.testdata.infiles['ps'], calmode='tsys',
                       outfile=self.testdata.tempfiles['applyfile'])
        if not args.get('outfile'):
            args['outfile'] = self.testdata.tempfiles['outfile']
        if whether is SUCCESS:
            self.assertIsNotNone(desired)
            self.positive_test(self.task, args, desired)
        else:
            self.negative_test(self.task, args)

    def test_sdcal_ps(self):
        """Test sdcal(calmode=ps)."""
        self.__test_sdcal(self.fail,
                          args={'infile': self.testdata.infiles['ps'],
                                'calmode': 'ps'})

    def test_sdcal_otfraster(self):
        """Test sdcal(calmode=otfraster)."""
        self.__test_sdcal(self.success,
                          args={'infile': self.testdata.infiles['otf'],
                                'calmode': 'otfraster'},
                          desired={'intent': 'OBSERVE_TARGET#ON_SOURCE'})

    def test_sdcal_otfraster_not_override(self):
        """Test sdcal(calmode=otfraster, intent='..')."""
        self.__test_sdcal(self.fail,
                          args={'infile': self.testdata.infiles['otf'],
                                'calmode': 'otfraster',
                                'intent': 'OBSERVE_TARGET#ON_SOURCE'})

    def test_sdcal_otf(self):
        """Test sdcal(calmode=otf)."""
        self.__test_sdcal(self.success,
                          args={'infile': self.testdata.infiles['otf'],
                                'calmode': 'otf'},
                          desired={'intent': 'OBSERVE_TARGET#ON_SOURCE'})

    def test_sdcal_otf_not_override(self):
        """Test sdcal(calmode=otf, intent='..')."""
        self.__test_sdcal(self.fail,
                          args={'infile': self.testdata.infiles['otf'],
                                'calmode': 'otf',
                                'intent': 'OBSERVE_TARGET#ON_SOURCE'})

    def test_sdcal_apply(self):
        """Test sdcal(calmode=apply)."""
        self.__test_sdcal(self.success,
                          args={'infile': self.testdata.infiles['ps'],
                                'calmode': 'apply',
                                'applytable': self.testdata.tempfiles['applyfile']},
                          desired={'interp': ''},
                          prepare_data=True)

    def test_sdcal_apply_not_override(self):
        """Test sdcal(calmode=apply, interp='..')."""
        self.__test_sdcal(self.fail,
                          args={'infile': self.testdata.infiles['ps'],
                                'calmode': 'apply',
                                'applytable': self.testdata.tempfiles['applyfile'],
                                'interp': 'nearest'},
                          prepare_data=True)

    def test_sdcal_ps_apply(self):
        """Test sdcal(calmode=ps,apply)."""
        self.__test_sdcal(self.success,
                          args={'infile': self.testdata.infiles['ps'],
                                'calmode': 'ps,apply'},
                          desired={'applytable': '',
                                   'interp': ''})

    def test_sdcal_ps_apply_not_override(self):
        """Test sdcal(calmode=ps,apply) with applytabe and interp are not null."""
        self.__test_sdcal(self.fail,
                          args={'infile': self.testdata.infiles['ps'],
                                'calmode': 'ps,apply',
                                'applytable': self.testdata.tempfiles['applyfile'],
                                'interp': 'nearest'},
                          prepare_data=True)

    def test_sdcal_tsys_apply(self):
        """Test sdcal(calmode=tsys,apply)."""
        self.__test_sdcal(self.success,
                          args={'infile': self.testdata.infiles['ps'],
                                'calmode': 'tsys,apply'},
                          desired={'applytable': '',
                                   'interp': ''})

    def test_sdcal_tsys_apply_not_override(self):
        """Test sdcal(calmode=tsys,apply) with applytabe and interp are not null."""
        self.__test_sdcal(self.fail,
                          args={'infile': self.testdata.infiles['ps'],
                                'calmode': 'tsys,apply',
                                'applytable': self.testdata.tempfiles['applyfile'],
                                'interp': 'nearest'},
                          prepare_data=True)

    def test_sdcal_ps_tsys_apply(self):
        """Test sdcal(calmode=ps,tsys,apply)."""
        self.__test_sdcal(self.success,
                          args={'infile': self.testdata.infiles['ps'],
                                'calmode': 'ps,tsys,apply'},
                          desired={'applytable': '',
                                   'interp': ''})

    def test_sdcal_ps_tsys_apply_not_override(self):
        """Test sdcal(calmode=ps,tsys,apply) with applytabe and interp are not null."""
        self.__test_sdcal(self.fail,
                          args={'infile': self.testdata.infiles['ps'],
                                'calmode': 'ps,tsys,apply',
                                'applytable': self.testdata.tempfiles['applyfile'],
                                'interp': 'nearest'},
                          prepare_data=True)

    def test_sdcal_otfraster_apply(self):
        """Test sdcal(calmode=otfraster,apply)."""
        self.__test_sdcal(self.success,
                          args={'infile': self.testdata.infiles['otf'],
                                'calmode': 'otfraster,apply'},
                          desired={'applytable': '',
                                   'interp': '',
                                   'intent': 'OBSERVE_TARGET#ON_SOURCE'})

    def test_sdcal_otfraster_apply_not_override(self):
        """Test sdcal(calmode=otfraster,apply) with applytabe/interp/intent are not null."""
        self.__test_sdcal(self.fail,
                          args={'infile': self.testdata.infiles['otf'],
                                'calmode': 'otfraster,apply',
                                'applytable': self.testdata.tempfiles['applyfile'],
                                'interp': 'nearest',
                                'intent': 'OBSERVE_TARGET#ON_SOURCE'},
                          prepare_data=True)

    def test_sdcal_otfraster_tsys_apply(self):
        """Test sdcal(calmode=otfraster,tsys,apply)."""
        self.__test_sdcal(self.success,
                          args={'infile': self.testdata.infiles['otf'],
                                'calmode': 'otfraster,tsys,apply'},
                          desired={'applytable': '',
                                   'interp': ''})

    def test_sdcal_otfraster_tsys_apply_not_override(self):
        """Test sdcal(calmode=otfraster,tsys,apply) with applytabe and interp are not null."""
        self.__test_sdcal(self.fail,
                          args={'infile': self.testdata.infiles['otf'],
                                'calmode': 'otfraster,tsys,apply',
                                'applytable': self.testdata.tempfiles['applyfile'],
                                'interp': 'nearest'},
                          prepare_data=True)

    def test_sdcal_otf_tsys_apply(self):
        """Test sdcal(calmode=otf,tsys,apply)."""
        self.__test_sdcal(self.success,
                          args={'infile': self.testdata.infiles['otf'],
                                'calmode': 'otf,tsys,apply'},
                          desired={'applytable': '',
                                   'interp': ''})

    def test_sdcal_otf_tsys_apply_not_override(self):
        """Test sdcal(calmode=otf,tsys,apply) with applytabe and interp are not null."""
        self.__test_sdcal(self.fail,
                          args={'infile': self.testdata.infiles['otf'],
                                'calmode': 'otf,tsys,apply',
                                'applytable': self.testdata.tempfiles['applyfile'],
                                'interp': 'nearest'},
                          prepare_data=True)


class SDFitTest(CasaxmlutilTest):
    """Test the constraints of sdfit."""

    testdata = Testdata(ctsys.resolve('unittest/sdfit'),
                        {'timebin': 'sdfit_tave.ms'},
                        {})

    task = sdfit
    success = SUCCESS
    fail = FAIL

    def __test_sdfit(self, whether: bool, args: dict, desired: dict=None):
        """Test sdfit with arguments and desired output."""
        if whether is SUCCESS:
            self.positive_test(self.task, args, desired)
        else:
            self.negative_test(self.task, args)

    def test_sdfit_timebin(self):
        """Test sdfit(timebin='1s'). If timebin is specified a value, then timespan is overridden by ''."""
        self.__test_sdfit(self.success,
                          args={'infile': self.testdata.infiles['timebin'],
                                'datacolumn': 'float_data',
                                'nfit': [1], 'pol': 'XX',
                                'timebin': '1s'},
                          desired={'timespan': ''})

    def test_sdfit_timebin_is_empty(self):
        """Test sdfit(timebin='')."""
        self.__test_sdfit(self.fail,
                          args={'infile': self.testdata.infiles['timebin'],
                                'datacolumn': 'float_data',
                                'nfit': [1], 'pol': 'XX',
                                'timebin': ''})


import importlib
casashell = importlib.find_loader('casashell')
if casashell is not None:
    from casashell.private.sdcal import sdcal as s_sdcal
    from casashell.private.sdfit import sdfit as s_sdfit

    class CasashellSDCalTest(SDCalTest):
        """Test that all constraint applying processes of casatasks.private.sdcal are not executed."""

        task = s_sdcal
        success = FAIL

    class CasashellSDFitTest(SDFitTest):
        """Test that all constraint applying processes of casatasks.private.sdfit are not executed."""

        task = s_sdfit
        success = FAIL

else:
    casalog.post('could not load tasks in casashell.private')


if __name__ == '__main__':
    unittest.main()
