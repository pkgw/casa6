import unittest
import os
import json
import shutil

from casatools import ctsys 
from casatestutils.stakeholder import almastktestutils
from casatestutils import compare

class Test_almastkutils(unittest.TestCase):
    
    @classmethod
    def setUpClass(cls):
        ''' common input data '''
        cls.datapath = ctsys.resolve('stakeholder/alma')
        cls.oldtestscr = 'test_stk_alma_pipeline_imaging_old.py'
        cls.curjson = 'test_stk_alma_pipeline_imaging_exp_dicts.json'
        cls.stdcubejson = 'test_standard_cube_cur_stats.json'
        cls.updatedjson = 'update_exp_subdicts.json'
       
    def setUp(self):
        shutil.copy(os.path.join(self.datapath,self.stdcubejson),self.stdcubejson)
        shutil.copy(os.path.join(self.datapath,self.curjson),self.curjson)
        shutil.copy(os.path.join(self.datapath,self.updatedjson), self.updatedjson)

    def compareDict(self, indict, refdict):
        ''' compare metric dictionaries, allow the comparison on vlaues between exp_dicts with extra flags and current metric dicts '''
        nfailcontent=0
        for key in indict:
            if key in refdict:
                if isinstance(indict[key],dict):
                    # additional dict level
                    for mt in indict[key]:
                        if mt in refdict[key]:
                            if isinstance(refdict[key][mt],list) and not isinstance(indict[key][mt],list):
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
                   if isinstance(refdict[key],list) and not isinstance(indict[key],list):
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

        compres = compare.compare_dictionaries(retdict,refdict['test_mosaic_cube'])
        self.assertTrue(compres)

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
        almastktestutils.update_expdict_subset(self.curjson, self.updatedjson, 'CAS-999999')
        outjsonfile = os.path.splitext(self.curjson)[0]+ '_update.json'
        print(outjsonfile)
        with open(outjsonfile) as f:
             expdict = json.load(f)
        
        # updatedjson contains updates for test_standard_cube and test_standard_mfs (sumwt values)
        self.assertEqual(expdict['test_standard_cube']['exp_sumwt_stats']['max_val'][1],18.0)
        self.assertEqual(expdict['test_standard_cube']['exp_sumwt_stats']['min_val'][1],17.0)
        self.assertEqual(expdict['test_standard_cube']['exp_sumwt_stats']['im_rms'][1],17.5)
        self.assertEqual(expdict['test_standard_mfs']['exp_sumwt_stats']['max_val'][1],6000400)
        self.assertEqual(expdict['test_standard_mfs']['exp_sumwt_stats']['min_val'][1],6000400)
        self.assertEqual(expdict['test_standard_mfs']['exp_sumwt_stats']['im_rms'][1],6000400)


    def test_compare_expdictjson(self):
        ''' Test a fuction to compare the two fiducial value json files '''

        # make another copy of a local copy of the fiducial value json 
        copyofthejson = os.path.splitext(self.curjson)[0]+'_copy.json'
        modcopyofthejson = os.path.splitext(self.curjson)[0]+'_copy.mod.json'
        shutil.copy(self.curjson, copyofthejson)

        # case 1: identical
        ret=almastktestutils.compare_expdictjson(self.curjson,copyofthejson)
        msg = 'The two json files are identical'

        self.assertTrue(msg in ret)

        # case 2: one of metric value is different
        with open(copyofthejson, 'r') as f:
            expdict2 = json.load(f)

        expdict2['test_standard_cube']['exp_im_stats']['im_rms'][1]=0.250
        with open(modcopyofthejson, 'w') as f2:
            json.dump(expdict2,f2)

        ret2=almastktestutils.compare_expdictjson(self.curjson,modcopyofthejson)
        refdict = {'test_standard_cube': {'exp_im_stats': {'im_rms': 
                                                               {'msg': 'diff in value(s)', 
                                                                'json1': [False, 0.143986095161], 
                                                                'json2': [False, 0.25]}}}}
        self.assertTrue(compare.compare_dictionaries(ret2,refdict))
