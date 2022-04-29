import unittest
import os
import json

from casatools import ctsys 
from casatestutils.stakeholder import almastktestutils

class TestAlmastkUtils(unittest.TestCase):
    def compareDict(self, indict, refdict):
        nfailcontent=0
        for key in indict:
            if key in refdict:
                if indict[key] == refdict[key][1]:
                    pass
                else:
                    nfailcontent += 1
            else:
                nfailcontent += 1
        return nfailcontent
                    
    def test_extract_eOxpdict(self):
        srcpath = ""
        almastktestutils.extract_expdict(testsrcpath=srcpath)

        self.assertTrue(os.path.exists('test_standard_cube_exp_dicts.json'))
        self.assertTrue(os.path.exists('test_standard_cube_eph_exp_dicts.json'))
        self.assertTrue(os.path.exists('test_standard_mfs_exp_dicts.json'))
        self.assertTrue(os.path.exists('test_standard_mfs_eph_exp_dicts.json'))
        self.assertTrue(os.path.exists('test_standard_mtmfs_exp_dicts.json'))
        self.assertTrue(os.path.exists('test_standard_mtmfs_eph_exp_dicts.json'))
        self.assertTrue(os.path.exists('test_standard_cal_exp_dicts.json'))
        self.assertTrue(os.path.exists('test_standard_cal_eph_exp_dicts.json'))
        self.assertTrue(os.path.exists('test_mosaic_cube_exp_dicts.json'))
        self.assertTrue(os.path.exists('test_mosaic_cube_eph_exp_dicts.json'))
        self.assertTrue(os.path.exists('test_mosaic_mfs_exp_dicts.json'))
        self.assertTrue(os.path.exists('test_mosaic_mfs_eph_exp_dicts.json'))
        self.assertTrue(os.path.exists('test_mosaic_mtmfs_exp_dicts.json'))
        self.assertTrue(os.path.exists('test_mosaic_mtmfs_eph_exp_dicts.json'))

    def test_create_expdict_jsonfile(self):
        injsonfile = 'test_standard_cube_tst.json'
        template = 'test_standard_cube_exp_dicts.json'
        outfile = 'test_standard_cube_new_exp_dicts.json'
        almastktestutils.create_expdict_jsonfile(injsonfile, template, outfile)

    def test_update_expdict_jsonfile(self):
        newexpdictlist = ['../data/std_cube_cur_expdicts.json']
        expdictjsonfile = '../data/test_stk_alma_pipeline_imaging_exp_dicts.json'
        almastktestutils.update_expdict_jsonfile(newexpdictlist, expdictjsonfile)

        contentcheck=False
        keycheck=False
        with open('../data/test_stk_alma_pipeline_imaging_exp_dicts_update.json') as f:
            alldict = json.load(f)
            with open(newexpdictlist[0]) as g:
                refdict = json.load(g)
            if 'test_standard_cube' in alldict:
                keycheck = True
                compres = self.compareDict(alldict['test_standard_cube'], refdict)

        self.assertTrue(keycheck)
        self.assertTrue(compres)


