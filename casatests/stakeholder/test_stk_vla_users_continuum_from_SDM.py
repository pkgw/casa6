##########################################################################
##########################################################################
# test_stk_vla_users_continuum_from_SDM.py
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
# https://open-jira.nrao.edu/browse/CASR-662
# https://open-jira.nrao.edu/browse/CAS-14036
#
#
##########################################################################

'''
VLA Users Stakeholder Continuum Test

Based on 3C 391 CASAguide Continuum Tutorial for CASA Version 6.2.0
https://casaguides.nrao.edu/index.php?title=VLA_Continuum_Tutorial_3C391-CASA6.2.0


Test list
- flag fraction and MS nrows after import, split and editing steps
- antenna position corrections retrieved by gencal
- fluxscale input table and bootstrapping results
- imstat values from the final pbcor image
- imstat values from the final residual image

'''

##########################################################################
##########################################################################

# Imports #
import os
import unittest
import numpy as np
import shutil
from glob import glob

from casatestutils import generate_weblog
from casatestutils import add_to_dict

from casatools import ctsys
from casatools import table as tbtool
from casatasks import casalog, importasdm, listobs, flagdata, gencal, setjy, gaincal, bandpass, fluxscale, applycal, split, statwt, tclean, impbcor, imstat
from casatasks.private.parallel.parallel_task_helper import ParallelTaskHelper

ctsys_resolve = ctsys.resolve

# location of data
data_path = ctsys_resolve('stakeholder/vla/')

clean_mask = "#CRTFv0 CASA Region Text Format version 0\npoly [[18:49:21.55474, -000.52.32.9190], [18:49:29.15960, -000.53.10.9387], [18:49:31.21928, -000.53.48.9583], [18:49:31.53618, -000.54.38.8592], [18:49:35.33866, -000.55.02.6210], [18:49:40.09177, -000.55.21.6298], [18:49:40.72560, -000.56.25.7880], [18:49:41.04265, -000.58.36.4810], [18:49:38.82455, -000.59.14.5014], [18:49:31.37791, -000.59.40.6415], [18:49:23.29749, -000.59.21.6320], [18:49:24.08969, -000.58.17.4736], [18:49:23.29750, -000.57.56.0874], [18:49:15.37558, -000.57.53.7105], [18:49:11.73154, -000.57.13.3138], [18:49:11.41470, -000.56.37.6701], [18:49:16.16792, -000.54.07.9680], [18:49:19.65353, -000.52.42.4239]] coord=J2000, corr=[I]"


##############################################
test_dict = {}
####    Tests     ####
class Test_vla_users_continuum(unittest.TestCase):

    def setUp(self):
        self._mytb = tbtool()
        self.refversion='6.2.0.124'
        self.sdmfile = 'TDEM0001_sb1218006_1.55310.33439732639'
        os.symlink(data_path+self.sdmfile, self.sdmfile)
        self.maskfile = '3c391_clean_mask.crtf'
        # os.symlink(data_path+self.maskfile, self.maskfile)
        self.writeMask()
        if ParallelTaskHelper.isMPIEnabled():
            self.parallel = True

    def tearDown(self):
        generate_weblog("test_stk_vla_users_continuum_from_SDM",test_dict)
        self._mytb.done()
        self.delData()

    def writeMask(self):
        with open(self.maskfile,'w') as out1:
            out1.write(clean_mask)

    def delData(self):
        os.unlink( self.sdmfile )
        os.unlink( self.maskfile )
        os.system('rm -rf '+self.sdmfile+'.ms*')
        del_files = glob('3c391_ctm*')
        for f in del_files:
            shutil.rmtree(f)

    def test_vla_users_continuum_fromSDM(self):
        """VLA Stakeholders tsts: Based on 3C 391 CASAguide Continuum Tutorial"""
        sdmname = self.sdmfile
        msname = '3c391_ctm_mosaic_10s_spw0.ms'
        msname_split = '3c391_ctm_mosaic_spw0.ms'

        test_name = self._testMethodName
        report=[]

        ## Data Import
        importasdm(asdm=sdmname,
                vis=sdmname+'.ms', createmms=False,
                ocorr_mode='co', lazy=False, asis='Receiver CalAtmosphere',
                process_caldevice=True, process_pointing=True, savecmds=True,
                outfile=sdmname+'.flagonline.txt',
                overwrite=False, bdfflags=False, with_pointing_correction=True)

        add_to_dict(self, output = test_dict, dataset = sdmname+'.ms')

        flagdata(vis=sdmname+'.ms', mode='list',
                inpfile=sdmname+'.flagonline.txt',
                tbuff=0.0, action='apply', flagbackup=False, savepars=False)

        flagdata(vis=sdmname+'.ms', mode='clip', clipzeros=True, 
                action='apply', flagbackup=False, savepars=False)

        flagdata(vis=sdmname+'.ms', mode='shadow', tolerance=0.0,
                action='apply', flagbackup=False, savepars=False)

        split(vis=sdmname+'.ms', outputvis=msname, keepflags=False, 
                spw='0', timebin='10s', datacolumn='DATA')

        add_to_dict(self, output = test_dict, dataset = msname)


        ## Data Editing
        flagdata(vis=msname, flagbackup=True, mode='manual', scan='1')

        flagdata(vis=msname, flagbackup=True, mode='manual', antenna='ea13,ea15')
        
        flagdata(vis=msname, mode='quack', quackinterval=10.0, quackmode='beg')

        flagsummary = flagdata(vis=msname, mode='summary')           
        result, expected, tolerance = 100* flagsummary['spw']['0']['flagged']/flagsummary['spw']['0']['total'], 20.558870499248552, 1.0
        msg = f"Expected flag percentage of {expected}, got {result}, tolerance {tolerance}"
        print(msg); report.append(msg)
        self.assertAlmostEqual(result, expected, delta=tolerance , msg=msg )


        listobs_output = listobs(vis=msname)
        result, expected = listobs_output['numrecords'], 850358 #845379
        msg = f"Expected {expected} rows in split MS, got {result}, tolerance exact"
        print(msg); report.append(msg)
        self.assertTrue(result==expected, msg=msg )


        ## Calibration
        gencal(vis=msname,caltable='3c391_ctm_mosaic_10s_spw0.antpos',caltype='antpos')

        # TEST: compare contents of gencal table
        self._mytb.open('3c391_ctm_mosaic_10s_spw0.antpos') 
        ant_array=self._mytb.getcol('FPARAM').squeeze().T 
        self._mytb.close()
        nonzero_ant_array=ant_array[np.logical_not(np.all( ant_array==0, axis=1 ))]             

        result, expected = nonzero_ant_array.shape[0], 14
        msg = f"Expected {expected} nonzero antenna position corrections, got {result}, tolerance exact"
        print(msg); report.append(msg)
        self.assertTrue(result==expected, msg=msg )

        result, expected, tolerance = np.min(nonzero_ant_array), -0.0257, 0.00001
        msg = f"Expected minimum antenna position correction of {expected}, got {result}, tolerance {tolerance}"
        print(msg); report.append(msg)
        self.assertAlmostEqual(result, expected, delta=tolerance , msg=msg )

        result, expected, tolerance = np.max(nonzero_ant_array), 0.0045, 0.00001
        msg = f"Expected maximum antenna position correction of {expected}, got {result}, tolerance {tolerance}"
        print(msg); report.append(msg)
        self.assertAlmostEqual(result, expected, delta=tolerance , msg=msg )


        setjy(vis=msname,field='J1331+3030',standard='Perley-Butler 2017',
                model='3C286_C.im',usescratch=True,scalebychan=True,spw='')

        flagdata(vis=msname,
                flagbackup=True, mode='manual', antenna='ea05')

        gaincal(vis=msname, caltable='3c391_ctm_mosaic_10s_spw0.G0', 
                field='J1331+3030', refant='ea21', spw='0:27~36', calmode='p', solint='int', 
                minsnr=5, gaintable=['3c391_ctm_mosaic_10s_spw0.antpos'])

        gaincal(vis=msname,caltable='3c391_ctm_mosaic_10s_spw0.K0', 
                field='J1331+3030',refant='ea21',spw='0:5~58',gaintype='K', 
                solint='inf',combine='scan',minsnr=5,
                gaintable=['3c391_ctm_mosaic_10s_spw0.antpos',
                        '3c391_ctm_mosaic_10s_spw0.G0'])

        bandpass(vis=msname,caltable='3c391_ctm_mosaic_10s_spw0.B0',
                field='J1331+3030',spw='',refant='ea21',combine='scan', 
                solint='inf',bandtype='B',
                gaintable=['3c391_ctm_mosaic_10s_spw0.antpos',
                            '3c391_ctm_mosaic_10s_spw0.G0',
                            '3c391_ctm_mosaic_10s_spw0.K0'])

        gaincal(vis=msname,caltable='3c391_ctm_mosaic_10s_spw0.G1',
                field='J1331+3030',spw='0:5~58',
                solint='inf',refant='ea21',gaintype='G',calmode='ap',solnorm=False,
                gaintable=['3c391_ctm_mosaic_10s_spw0.antpos',
                        '3c391_ctm_mosaic_10s_spw0.K0',
                        '3c391_ctm_mosaic_10s_spw0.B0'],
                interp=['','','nearest'])

        gaincal(vis=msname,caltable='3c391_ctm_mosaic_10s_spw0.G1',
                field='J1822-0938',
                spw='0:5~58',solint='inf',refant='ea21',gaintype='G',calmode='ap',
                gaintable=['3c391_ctm_mosaic_10s_spw0.antpos',
                        '3c391_ctm_mosaic_10s_spw0.K0',
                        '3c391_ctm_mosaic_10s_spw0.B0'],
                append=True)

        myscale = fluxscale(vis=msname,
                caltable='3c391_ctm_mosaic_10s_spw0.G1', 
                fluxtable='3c391_ctm_mosaic_10s_spw0.fluxscale1', 
                reference='J1331+3030',
                transfer=['J1822-0938'],
                incremental=False)


        # TEST: compare specific values in return dictionary
        result, expected = myscale['1']['0']['numSol'][0], 46
        msg = f"Expected {expected} input fluxboot solutions, got {result}, tolerance exact"
        print(msg); report.append(msg)
        self.assertTrue(result==expected, msg=msg )

        result, expected, epsilon = myscale['1']['0']['fluxd'][0], 2.2973563472684653, 0.04
        msg = f"Expected J1822-0938 bootstrapped flux of {expected}, got {result}, tolerance {100*epsilon}%"
        print(msg); report.append(msg)
        self.assertAlmostEqual(result, expected, delta=expected*epsilon, msg=msg )


        applycal(vis=msname,
                field='J1331+3030',
                gaintable=['3c391_ctm_mosaic_10s_spw0.antpos', 
                            '3c391_ctm_mosaic_10s_spw0.fluxscale1',
                            '3c391_ctm_mosaic_10s_spw0.K0',
                            '3c391_ctm_mosaic_10s_spw0.B0'],
                gainfield=['','J1331+3030','',''], 
                interp=['','nearest','',''],
                calwt=False)

        applycal(vis=msname,
                field='J1822-0938',
                gaintable=['3c391_ctm_mosaic_10s_spw0.antpos', 
                            '3c391_ctm_mosaic_10s_spw0.fluxscale1',
                            '3c391_ctm_mosaic_10s_spw0.K0',
                            '3c391_ctm_mosaic_10s_spw0.B0'],
                gainfield=['','J1822-0938','',''], 
                interp=['','nearest','',''],
                calwt=False)

        applycal(vis=msname,
                field='2~8',
                gaintable=['3c391_ctm_mosaic_10s_spw0.antpos', 
                            '3c391_ctm_mosaic_10s_spw0.fluxscale1',
                            '3c391_ctm_mosaic_10s_spw0.K0',
                            '3c391_ctm_mosaic_10s_spw0.B0'],
                gainfield=['','J1822-0938','',''], 
                interp=['','linear','',''],
                calwt=False)


        ## Imaging
        split(vis=msname,outputvis=msname_split,
                datacolumn='corrected',field='2~8',correlation='RR,LL')

        add_to_dict(self, output = test_dict, dataset = msname_split)

        statwt(vis=msname_split,datacolumn='data')

        tclean(vis=msname_split,imagename='3c391_ctm_spw0_multiscale',
                field='',spw='',
                specmode='mfs',
                niter=20000,
                gain=0.1, threshold='1.0mJy',
                gridder='mosaic',
                deconvolver='multiscale',
                scales=[0, 5, 15, 45], smallscalebias=0.9,
                interactive=False,
                imsize=[480,480], cell=['2.5arcsec','2.5arcsec'],
                stokes='I',
                weighting='briggs',robust=0.5,
                pbcor=False,
                mask=self.maskfile,
                savemodel='none',
               parallel=self.parallel)

        impbcor(imagename='3c391_ctm_spw0_multiscale.image',pbimage='3c391_ctm_spw0_multiscale.pb',
                outfile='3c391_ctm_spw0_multiscale.pbcorimage')

        mystat = imstat(imagename='3c391_ctm_spw0_multiscale.pbcorimage')

        # TESTS: compare values in return dictionary

        # location of image peak
        result, expected = mystat['maxpos'][:2], np.array([288,256])
        msg = f"Expected image peak at pixels {expected}, got {result}, tolerance exact"
        print(msg); report.append(msg)
        self.assertTrue( np.all(result==expected), msg=msg )

        # value of image peak
        result, expected, epsilon = mystat['max'][0], 0.15553903579711914, 0.04
        msg = f"Expected pbcor image peak of {expected}, got {result}, tolerance {100*epsilon}%"
        print(msg); report.append(msg)
        self.assertAlmostEqual(result, expected, delta=expected*epsilon, msg=msg )


        mystat = imstat(imagename='3c391_ctm_spw0_multiscale.residual', region=self.maskfile)

        # TESTS: compare values in return dictionary

        # max residual inside mask
        result, expected, epsilon = mystat['max'][0], 0.0010077510960400105, 0.10
        msg = f"Expected peak residual of {expected}, got {result}, tolerance {100*epsilon}%"
        print(msg); report.append(msg)
        self.assertAlmostEqual(result, expected, delta=expected*epsilon, msg=msg )

        # peak residual inside mask
        result, expected, epsilon = mystat['rms'][0], 0.000542584067811474, 0.10
        msg = f"Expected rms residual of {expected}, got {result}, tolerance {100*epsilon}%"
        print(msg); report.append(msg)
        self.assertAlmostEqual(result, expected, delta=expected*epsilon, msg=msg )


        test_dict[test_name]['report'] = '\n'.join(report)
        test_dict[test_name]['images'] = []

# main
if __name__ == '__main__':
    unittest.main()
