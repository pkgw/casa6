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
                if type(indict[key]==dict):
                    # additional dict level
                    for mt in indict[key]:
                        if mt in refdict[key]:
                            if type(refdict[key][mt])==list and type(indict[key][mt])!=list:
                                if indict[key][mt] == refdict[key][mt][1]:
                                    pass
                                else:
                                    nfailcontent += 1
                            else:
                                if indict[key][mt] == refdict[key][mt]:
                                    pass
                                else:
                                    nfailcontent += 1
                        else:
                            nfailcontent += 1
                else: 
                   if type(refdict[key])==list and type(indict[key])!=list:
                      if indict[key] == refdict[key][1]:
                          pass
                      else:
                          nfailcontent += 1
                   else:
                      if indict[key]==refdict[key]:
                          pass 
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

    def test_read_testcase_expdicts(self):
        ''' Test a fuction to read the combined metric value json and returns a dictionary of the values for a testcase'''
        # get current version info from the fiducial value json
        with open (os.path.join(self.datapath, self.curjson)) as f:
             refdict = json.load(f)
        casaversion = ''
        if 'casa_version' in refdict:
           casaversion = refdict['casa_version']
           self.assertNotEqual(casaversion,'')
 
        retdict = almastktestutils.read_testcase_expdicts(os.path.join(self.datapath,self.curjson), 'test_mosaic_cube', casaversion)

        compres = self.compareDict(retdict,refdict['test_mosaic_cube'])
        self.assertTrue(compres==0)

    def test_read_testcase_expdicts_mismatchCasaVersion(self):
        ''' Test a fuction to read stored fiducial value json: CASA version check '''
        try:
            almastktestutils.read_testcase_expdicts(os.path.join(self.datapath,self.curjson), 'test_mosaic_cube','6.2.0')
        except Exception as e:
            msg = 'Mismatch in the fiducial data file version. The testcase expects fiducial values based on the CASA 6.2.0'
            self.assertTrue(msg in str(e))


    def test_create_expdict_jsonfile(self):
        outfile = 'test_standard_cube_new_exp_dicts.json'
        almastktestutils.create_expdict_jsonfile(self.stdcubejson, os.path.join(self.datapath,self.curjson), outfile)

        has_test_standard_cube = False
        nmismatch = 0
        outputexist = os.path.exists(outfile)
        self.assertTrue(outputexist)
        if outputexist:
            with open(outfile) as f, open(os.path.join(self.datapath,self.curjson)) as g:
                outdict = json.load(f)
                fulldict = json.load(g)

            if 'test_standard_cube' in outdict:
                has_test_standard_cube = True
                for mt in outdict['test_standard_cube']:
                    if mt in fulldict['test_standard_cube']:
                       for mc in outdict['test_standard_cube'][mt]:
                           if mc in fulldict['test_standard_cube']:
                               if outdict['test_stanard_cube'][mt][mc] == fulldict['test_standard_cube'][mt][mc]:
                                   pass
                               else: 
                                   nmismatch += 1                
 
            self.assertTrue(has_test_standard_cube)
            self.assertEqual(nmismatch,0)

  
    def test_update_expdict_jsonfile(self):
        ''' Test a fuction that updates the combined exp dict json using list of metric value json files of individual testcases '''
        newexpdictlist = [self.stdcubejson]
        expdictjsonfile = self.curjson 
        almastktestutils.update_expdict_jsonfile(newexpdictlist, expdictjsonfile)

        compres=False
        keycheck=False
        with open(expdictjsonfile.split('.json')[0]+'_update.json') as f:
            alldict = json.load(f)
            with open(newexpdictlist[0]) as g:
                refdict = json.load(g)
            if 'test_standard_cube' in alldict:
                keycheck = True
                compres = self.compareDict(alldict['test_standard_cube'], refdict)
            else:
                print("skipping content test")

        self.assertTrue(keycheck)
        self.assertTrue(compres==0)


    def test_update_expdict_subset(self):
        ''' Test a fuction to update only subset of metric values to the combined expdict json '''
         # update_expdict_subset(expjsonfile, newvaldictjson, jiranoforcomment='')
        pass


    def test_compare_expdictjson(self):
        pass
 
