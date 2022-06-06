import unittest
import os
import json
import shutil

from casatools import ctsys 
from casatestutils.stakeholder import almastktestutils

class Test_almastkutils(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        ''' common input data '''
        cls.datapath = ctsys.resolve('stakeholder/alma')
        cls.oldtestscr = 'test_stk_alma_pipeline_imaging_old.py'
        cls.curjson = 'test_stk_alma_pipeline_imaging_exp_dicts.json'
        cls.stdcubejson = 'test_standard_cube_cur_stats.json'
       
    def setUp(self):
        shutil.copy(os.path.join(self.datapath,self.stdcubejson),self.stdcubejson)
        shutil.copy(os.path.join(self.datapath,self.curjson),self.curjson)

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
                    
    def test_extract_subexpdict(self):
        ''' Test a fuction to extract a subset of exp dicts from a single testcase json file '''
        print("stdcubejon=",self.stdcubejson)
        metrickeys=['im_stats_dict','psf_stats_dict']
        # Use default for outjsonfile. The output json will be inputjsonfilename+'subDict.json'
        almastktestutils.extract_subexpdict(self.stdcubejson,keylist=metrickeys)
        outputname = self.stdcubejson.rstrip('.json') + '_subDict.json'
        print("output=",outputname)
        outputexist = os.path.exists(outputname)
        self.assertTrue(outputexist)
        if outputexist:
            with open(outputname) as f:
               subdict = json.load(f)
            self.assertTrue('test_standard_cube' in subdict)
            self.assertTrue('im_stats_dict' in subdict['test_standard_cube']) 

            self.assertTrue('psf_stats_dict' in subdict['test_standard_cube']) 

    def test_create_expdict_jsonfile(self):
        outfile = 'test_standard_cube_new_exp_dicts.json'
        almastktestutils.create_expdict_jsonfile(self.stdcubejson, os.path.join(self.datapath,self.curjson), outfile)

    def test_update_expdict_jsonfile(self):
        ''' Test a fuction that updates the combined exp dict json using list of metric value json files of individual testcases '''
        newexpdictlist = [self.stdcubejson]
        expdictjsonfile = self.curjson 
        almastktestutils.update_expdict_jsonfile(newexpdictlist, expdictjsonfile)

        contentcheck=False
        keycheck=False
        with open(expdictjsonfile.split('.json')[0]+'_update.json') as f:
            alldict = json.load(f)
            with open(newexpdictlist[0]) as g:
                refdict = json.load(g)
            if 'test_standard_cube' in alldict:
                keycheck = True
                compres = self.compareDict(alldict['test_standard_cube'], refdict)

        self.assertTrue(keycheck)
        self.assertTrue(compres)


    def test_update_expdict_subset(self):
        ''' Test a fuction to update only subset of metric values to the combined expdict json '''
         # update_expdict_subset(expjsonfile, newvaldictjson, jiranoforcomment='')
        pass 
