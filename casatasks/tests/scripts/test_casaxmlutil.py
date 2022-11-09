
import os
import shutil
import unittest

from casatasks.private.casaxmlutil import xml_constraints_injector
from casatasks import sdcal, casalog
from casatools import ctsys

casalogpath = casalog.logfile()
testdir = 'cxutest'

sdcal_datapath = ctsys.resolve('unittest/sdcal')
sdcal_infiles = {'ps': 'uid___A002_X6218fb_X264.ms.sel', 'otf': 'uid___A002_X6218fb_X264.ms.sel.otfraster'}
sdcal_outfile = 'sdcal.out'
sdcal_applyfile = 'apply.cal'


class casaxmlutil_test(unittest.TestCase):

    @xml_constraints_injector
    def dummy(self, *args, **kwargs):
        return args, kwargs

    @classmethod
    def setUpClass(cls):
        cls.curdir = os.getcwd()
        if os.path.exists(testdir):
            shutil.rmtree(testdir)
        os.mkdir(testdir)
        os.chdir(testdir)
        for v in sdcal_infiles.values():
            if os.path.exists(v):
                shutil.rmtree(v)
            shutil.copytree(os.path.join(sdcal_datapath, v), v)

    @classmethod
    def tearDownClass(cls):
        os.chdir(cls.curdir)
        # if os.path.exists(testdir):
        #     shutil.rmtree(testdir)

    def tearDown(self):
        casalog.setlogfile(casalogpath)
        if os.path.exists(sdcal_outfile):
            shutil.rmtree(sdcal_outfile)
        if os.path.exists(sdcal_applyfile):
            shutil.rmtree(sdcal_applyfile)

    def __load_logfile(self, logfile):
        with open(logfile, 'r') as fp:
            for line in fp:
                yield line.rstrip()

    def __check_log(self, logfile, msg):
        for line in self.__load_logfile(logfile):
            if msg in line:
                return True
        return False

    def test_sdcal_ps(self):
        logfile = 'sdcal_log_ps.txt'
        casalog.setlogfile(logfile)
        args = {'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                'calmode': 'ps'}
        sdcal(**args)
        self.assertFalse(self.__check_log(logfile, "overrode argument:"))

    def test_sdcal_otfraster(self):
        logfile = 'sdcal_log_otfraster.txt'
        casalog.setlogfile(logfile)
        args = {'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                'calmode': 'otfraster'}
        sdcal(**args)
        self.assertTrue(self.__check_log(logfile, "overrode argument: intent -> 'OBSERVE_TARGET#ON_SOURCE'"))

    def test_sdcal_otf(self):
        logfile = 'sdcal_log_otf.txt'
        casalog.setlogfile(logfile)
        args = {'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                'calmode': 'otf'}
        sdcal(**args)
        self.assertTrue(self.__check_log(logfile, "overrode argument: intent -> 'OBSERVE_TARGET#ON_SOURCE'"))

    def test_sdcal_otf_not_override(self):
        logfile = 'sdcal_log_otf_not_override.txt'
        casalog.setlogfile(logfile)
        args = {'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                'calmode': 'otf', 'intent': 'OBSERVE_TARGET#ON_SOURCE'}
        sdcal(**args)
        self.assertFalse(self.__check_log(logfile, "overrode argument:"))

    def test_sdcal_apply(self):
        sdcal(infile=sdcal_infiles['ps'], calmode='tsys', outfile=sdcal_applyfile)
        logfile = 'sdcal_log_apply.txt'
        casalog.setlogfile(logfile)
        args = {'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                'calmode': 'apply', 'applytable': sdcal_applyfile}
        sdcal(**args)
        self.assertTrue(self.__check_log(logfile, "overrode argument: interp -> ''"))

    def test_sdcal_ps_apply(self):
        logfile = 'sdcal_log_ps_apply.txt'
        casalog.setlogfile(logfile)
        args = {'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                'calmode': 'ps,apply'}
        sdcal(**args)
        self.assertTrue(self.__check_log(logfile, "overrode argument: applytable -> ''"))
        self.assertTrue(self.__check_log(logfile, "overrode argument: interp -> ''"))
        self.assertTrue(self.__check_log(logfile, "recursive task call"))

    def test_sdcal_tsys_apply(self):
        logfile = 'sdcal_log_tsys_apply.txt'
        casalog.setlogfile(logfile)
        args = {'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                'calmode': 'tsys,apply'}
        sdcal(**args)
        self.assertTrue(self.__check_log(logfile, "overrode argument: applytable -> ''"))
        self.assertTrue(self.__check_log(logfile, "overrode argument: interp -> ''"))
        self.assertTrue(self.__check_log(logfile, "recursive task call"))

    def test_sdcal_ps_tsys_apply(self):
        logfile = 'sdcal_log_ps_tsys_apply.txt'
        casalog.setlogfile(logfile)
        args = {'infile': sdcal_infiles['ps'], 'outfile': sdcal_outfile,
                'calmode': 'ps,tsys,apply'}
        sdcal(**args)
        self.assertTrue(self.__check_log(logfile, "overrode argument: applytable -> ''"))
        self.assertTrue(self.__check_log(logfile, "overrode argument: interp -> ''"))
        self.assertTrue(self.__check_log(logfile, "recursive task call"))

    def test_sdcal_otfraster_apply(self):
        logfile = 'sdcal_log_otfraster_apply.txt'
        casalog.setlogfile(logfile)
        args = {'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                'calmode': 'otfraster,apply'}
        sdcal(**args)
        self.assertTrue(self.__check_log(logfile, "overrode argument: applytable -> ''"))
        self.assertTrue(self.__check_log(logfile, "overrode argument: interp -> ''"))
        self.assertTrue(self.__check_log(logfile, "overrode argument: intent -> 'OBSERVE_TARGET#ON_SOURCE'"))
        self.assertTrue(self.__check_log(logfile, "recursive task call"))

    def test_sdcal_otfraster_tsys_apply(self):
        logfile = 'sdcal_log_otfraster_tsys_apply.txt'
        casalog.setlogfile(logfile)
        args = {'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                'calmode': 'otfraster,tsys,apply'}
        sdcal(**args)
        self.assertTrue(self.__check_log(logfile, "overrode argument: applytable -> ''"))
        self.assertTrue(self.__check_log(logfile, "overrode argument: interp -> ''"))
        self.assertTrue(self.__check_log(logfile, "recursive task call"))

    def test_sdcal_otf_apply(self):
        logfile = 'sdcal_log_otf_apply.txt'
        casalog.setlogfile(logfile)
        args = {'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                'calmode': 'otf,apply'}
        sdcal(**args)
        self.assertTrue(self.__check_log(logfile, "overrode argument: applytable -> ''"))
        self.assertTrue(self.__check_log(logfile, "overrode argument: interp -> ''"))
        self.assertTrue(self.__check_log(logfile, "overrode argument: intent -> 'OBSERVE_TARGET#ON_SOURCE'"))
        self.assertTrue(self.__check_log(logfile, "recursive task call"))

    def test_sdcal_otf_tsys_apply(self):
        logfile = 'sdcal_log_otf_tsys_apply.txt'
        casalog.setlogfile(logfile)
        args = {'infile': sdcal_infiles['otf'], 'outfile': sdcal_outfile,
                'calmode': 'otf,tsys,apply'}
        sdcal(**args)
        self.assertTrue(self.__check_log(logfile, "overrode argument: applytable -> ''"))
        self.assertTrue(self.__check_log(logfile, "overrode argument: interp -> ''"))
        self.assertTrue(self.__check_log(logfile, "recursive task call"))


if __name__ == '__main__':
    unittest.main()
