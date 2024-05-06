##########################################################################
# test_stk_srdp_restore_mode.py
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
# https://open-jira.nrao.edu/browse/CASR-517
# https://open-jira.nrao.edu/browse/CAS-14319
#
##########################################################################

'''
These tests will run importadsm to fill a CASA MS and restore the flag backup saved
in a previous run of the same project. The tests will check the flagging stats between running the 
calibration tasks to verify that the flags have not changed compared to the expected values 
from the previous run. The last step will create a calibrator image and check the flux density 
and beam values.

This script runs the restore mode of the pipeline to 3 use-cases: ALMA 12m, 7m, and VLA.
'''

import os
import shutil
import re
import unittest
from casatools import ctsys
from casatasks import importasdm, flagdata, flagmanager, applycal, tclean, imstat, imhead, hanningsmooth, statwt
from casatestutils import stats_dict
from casatestutils import generate_weblog
from casatestutils import add_to_dict
from casatestutils.imagerhelpers import TestHelpers
th = TestHelpers()

data_path = ctsys.resolve('stakeholder/srdp/')


class Test_srdp_base(unittest.TestCase):
    """ Base class for all tests """
    
    @classmethod
    def tearDownClass(cls) -> None:
       # Generate weblogs for the tests. An html file will be created in the local directory
        generate_weblog("srdp_restore_mode",test_dict)
        #return super().tearDownClass()
    
    def _extract_caltables(self, filename):
        # This pattern matches the 'caltable' key and extracts its value
        pattern = re.compile(r"caltable='(.*?)'")
        caltables = []

        with open(filename, 'r') as file:
            for line in file:
                match = pattern.search(line)
                if match:
                    # Extract the name of the caltable and add it to the list
                    caltables.append(match.group(1))
            unique_tables = set(caltables)

        return unique_tables

    def filter_report(self, report, showonlyfail=True):
        """ Function to filter the test report, the input report is expected to be a string with the newline code """
        ret = ''
        if showonlyfail:
            filter='Fail'
        else:
            filter='Pass'

        if report!='':
            testItems = report.split('\n')
            retitems=[]
            for testitem in testItems:
                if '[ check_ims ]' in testitem or '[ check_pixmask ]' in testitem or '[ check_val ]' in testitem:
                    if '( '+filter in testitem:
                        retitems.append(testitem)
            nfail = len(retitems)
            msg = str(nfail)+' individual test failure(s) '
            ret = '\n' + '\n'.join(retitems)
            ret += '\n' + msg
        return ret
        
##########   ALMA 12m test ############
test_dict = {}
class Test_srdp_alma_12m(Test_srdp_base):
    """ ALMA 12m dataset """

    def setUp(self) -> None:
        # Reference criteria using casa-6.5.4-9-pipeline-2023.1.0.125
        # flux density in Jy; major, minor in arcsec and position angle in deg
        self.exp_flags = {'flags_importasdm':3900108.0, 'flags_restored':159843772.0, 'flags_final':159843772.0}
        self.exp_total_vis = 1092709620.0
        self.exp_flux = 1.251
        self.exp_beam = {'major':0.543,'minor':0.345,"positionangle":-19.153}
        self.exp_tol = 0.01 # within 1%

        # Input data
        self.asdm_name = 'uid___A002_Xd3607d_X2b0d'
        self.restore_flags_tar = 'uid___A002_Xd3607d_X2b0d.ms.flagversions.tgz'
        # Reference callibrary has the names of the cal tables to apply
        self.reference_callibrary = 'uid___A002_Xd3607d_X2b0d.ms.s1.3.callibrary'

        # Output data
        self.output_ms = self.asdm_name + '.ms'
        self.restore_flags = 'uid___A002_Xd3607d_X2b0d.ms.flagversions'
        self.image_prefix = 'J0336+3218'

        # Get input data
        os.symlink(data_path + self.asdm_name, self.asdm_name)
        os.symlink(data_path + self.restore_flags_tar, self.restore_flags_tar)
        shutil.copy(data_path + self.reference_callibrary, self.reference_callibrary)
        self.caltables = self._extract_caltables(self.reference_callibrary)

        # Copy cal tables locally
        for ctable in self.caltables:
            shutil.copytree(data_path + ctable, ctable)

    def tearDown(self) -> None:
        os.unlink(self.asdm_name)
        os.unlink(self.restore_flags_tar)
        os.remove(self.reference_callibrary)
        for ctable in self.caltables:
            shutil.rmtree(ctable)
        
        os.system('rm -rf '+self.image_prefix + '*')
        shutil.rmtree(self.output_ms)
        shutil.rmtree(self.restore_flags)

    @stats_dict(test_dict)
    def test_alma_12m_restore(self):
        """ SRDP ALMA 12m pipeline restore mode: check that restored flags do not change in new MS """
        # Import the MS and do not apply the online flags
        # The BDF flags are saved with the flagbackup in the .flagversions file with name Original
        test_name = self._testMethodName
        report = []
        importasdm(asdm=self.asdm_name, vis=self.output_ms, createmms=False, ocorr_mode='ca', lazy=False,
                   asis='SBSummary ExecBlock Antenna Station Receiver Source CalAtmosphere CalWVR CalPointing',
                   process_caldevice=False, process_flags=True, applyflags=False, savecmds=False,
                   overwrite=False, bdfflags=True, with_pointing_correction=False)

        # Get the flags summary and compare with the expected values
        flags_importasdm = flagdata(vis=self.output_ms, mode='summary', name='before-flagmanager')

        # Check the total visibilities
        status, report0 = th.check_val(flags_importasdm['total'], self.exp_total_vis, exact=True, \
                               valname='total_visibilities')

        # Check the flags before flagmanager. There should be only the bdfflags in the MS
        status, report1 = th.check_val(flags_importasdm['flagged'], self.exp_flags['flags_importasdm'], exact=True, \
                               valname='flags_importasdm')

        # Restore the reference flags from the pipeline to the MS
        os.system('tar xf '+ self.restore_flags_tar)
        flagmanager(vis=self.output_ms, mode='restore', versionname='Pipeline_Final')

        # Get the flags summary and compare with expected values
        # At this point, the MS has the sum of flags from Original and Pipeline_Final
        flags_restored = flagdata(vis=self.output_ms, mode='summary', name='before-applycal')
        status, report2 = th.check_val(flags_restored['flagged'], self.exp_flags['flags_restored'], exact=True, \
                               valname='flags_restored')

        # Apply the calibration tables listed in the cal library
        applycal(vis=self.output_ms,
                 field='J0328+3139,J0336+3218,NGC1333_IRAS_4A,J0237+2848',
                 spw='19,21,23,25',
                 intent='CALIBRATE_BANDPASS#ON_SOURCE,OBSERVE_CHECK_SOURCE#ON_SOURCE,CALIBRATE_PHASE#ON_SOURCE,'\
                 'OBSERVE_TARGET#ON_SOURCE,CALIBRATE_FLUX#ON_SOURCE',
                 antenna='*&*', docallib=True,
                 callib=self.reference_callibrary,
                 applymode='calflagstrict', flagbackup=True)

        # Check that applycal did not flag anything
        flags_final = flagdata(vis=self.output_ms, mode='summary', name='applycal')
        status, report3 = th.check_val(flags_final['flagged'], self.exp_flags['flags_final'], exact=True, \
                               valname='flags_final')

        # Create the calibrator image and check the flux density and beam
        tclean(vis=self.output_ms,field='J0336+3218',imsize=256,cell='0.05arcsec',spw='19,21,23',
               imagename=self.image_prefix,niter=5000,nsigma=3.0)

        # Check calibrator peak flux density matches <1% of flux: 1.251 Jy 
        flux = imstat(imagename=self.image_prefix+'.image')['flux']
        status, report4 = th.check_val(flux, self.exp_flux, exact=False, epsilon=self.exp_tol, \
                               valname='flux')
 
        # Check that the beam maj/minor is within 1% of beam: 0.543" x 0.345" -19.153 deg
        beam = imhead(imagename=self.image_prefix+'.image')['restoringbeam']
        status, report5 = th.check_val(beam['major']['value'], self.exp_beam['major'], exact=False, epsilon=self.exp_tol, \
                                       valname='beam_major')

        status, report6 = th.check_val(beam['minor']['value'], self.exp_beam['minor'], exact=False, epsilon=self.exp_tol, \
                                       valname='beam_minor')
        
        status, report7 = th.check_val(beam['positionangle']['value'], self.exp_beam['positionangle'], exact=False, \
                                       epsilon=self.exp_tol, valname='positionangle') 
        
        report = report0 + report1 + report2 + report3 + report4 + report5 + report6 + report7
        failed = self.filter_report(str(report))
        test_dict[test_name]['report'] = report
        add_to_dict(self, output = test_dict, dataset = self.output_ms)

        self.assertTrue(th.check_final(pstr = report), msg = failed)

##########   ALMA 7m test ############
class Test_srdp_alma_7m(Test_srdp_base):
    """ ALMA 7m dataset """

    def setUp(self) -> None:
        # Reference criteria using casa-6.5.4-9-pipeline-2023.1.0.125
        # flagging before flagmanager matches (flagged and total); 'flagged': 2160.0, 'total': 378548640.0
        # flagging after flagmanager matches (flagged and total); 'flagged': 125159272.0, 'total': 378548640.0
        # flagging after applycal matches (flagged and total);  'flagged': 125159272.0,  'total': 378548640.0
        # calibrator peak flux density matches <1%, beam maj/minor within 1%; flux: 1.553 Jy; beam: 5.035" x 2.975" 78.692 deg

        self.exp_flags = {'flags_importasdm':2160.0, 'flags_restored':125159272.0, 'flags_applycal':125159272.0}
        self.exp_total_vis = 378548640.0
        self.exp_flux = 1.553 #Jy
        self.exp_beam = {'major':5.035,'minor':2.975,"positionangle":78.692}
        self.exp_tol = 0.01 # within 1%

        # Input data
        self.asdm_name = 'uid___A002_Xd341ff_X3a2f'
        self.restore_flags_tar = 'uid___A002_Xd341ff_X3a2f.ms.flagversions.tgz'
        # Reference callibrary has the names of the cal tables to apply
        self.reference_callibrary = 'uid___A002_Xd341ff_X3a2f.ms.s1.3.callibrary'

        # Output data
        self.output_ms = self.asdm_name + '.ms'
        self.restore_flags = 'uid___A002_Xd341ff_X3a2f.ms.flagversions'
        self.image_prefix = 'J0501-0159'

        # Get input data
        os.symlink(data_path + self.asdm_name, self.asdm_name)
        os.symlink(data_path + self.restore_flags_tar, self.restore_flags_tar)
        shutil.copy(data_path + self.reference_callibrary, self.reference_callibrary)
        self.caltables = self._extract_caltables(self.reference_callibrary)

        # Copy cal tables locally
        for ctable in self.caltables:
            shutil.copytree(data_path + ctable, ctable)

    def tearDown(self) -> None:
        os.unlink(self.asdm_name)
        os.unlink(self.restore_flags_tar)
        os.remove(self.reference_callibrary)
        for ctable in self.caltables:
            shutil.rmtree(ctable)
        
        os.system('rm -rf '+self.image_prefix + '*')
        shutil.rmtree(self.output_ms)
        shutil.rmtree(self.restore_flags)

    @stats_dict(test_dict)
    def test_alma_7m_restore(self):
        """ SRDP ALMA 7m pipeline restore mode: check that restored flags do not change in new MS """
        test_name = self._testMethodName
        report = []
        importasdm(asdm=self.asdm_name, vis=self.output_ms,createmms=False, ocorr_mode='ca', lazy=False, 
                   asis='SBSummary ExecBlock Antenna Station Receiver Source CalAtmosphere CalWVR CalPointing', 
                   process_caldevice=False, savecmds=False, overwrite=False, bdfflags=True, with_pointing_correction=False)

        # Check the flags summary and compare with expected values
        flags_importasdm = flagdata(vis=self.output_ms, mode='summary', name='before-flagmanager')

        # Check the total visibilities
        status, report0 = th.check_val(flags_importasdm['total'], self.exp_total_vis, exact=True, \
                               valname='total_visibilities')

        # Check the flags before flagmanager. There should be only the bdfflags in the MS
        status, report1 = th.check_val(flags_importasdm['flagged'], self.exp_flags['flags_importasdm'], exact=True, \
                               valname='flags_importasdm')
        
        # Restore the reference flags from the pipeline to the MS
        os.system('tar xf '+ self.restore_flags_tar)
        flagmanager(vis=self.output_ms, mode='restore', versionname='Pipeline_Final')

        # Get the flags summary and compare with expected values
        # At this point, the MS has the sum of flags from Original and Pipeline_Final
        flags_restored = flagdata(vis=self.output_ms, mode='summary', name='before-applycal')
        status, report2 = th.check_val(flags_restored['flagged'], self.exp_flags['flags_restored'], exact=True, \
                               valname='flags_restored')
        
       # Apply the calibration tables listed in the cal library
        applycal(vis=self.output_ms,
                field='HOPS-300,HOPS-325,HOPS-164,HOPS-166,HOPS-010,HOPS-185,HOPS-281,HOPS-050,J0501-0159,' \
                    'HOPS-404,HOPS-402,J0522-3627,HOPS-152,HOPS-337,HOPS-045,HOPS-042,HOPS-239,HOPS-397,' \
                    'HOPS-290,HOPS-247,HOPS-241,HOPS-240', spw='16,18,20,22',
                intent='CALIBRATE_BANDPASS#ON_SOURCE,OBSERVE_TARGET#ON_SOURCE,CALIBRATE_PHASE#ON_SOURCE,' \
                'CALIBRATE_FLUX#ON_SOURCE', antenna='*&*', docallib=True,
                callib=self.reference_callibrary, applymode='calflagstrict', flagbackup=True)
        
        # Check the flags after applycal
        flags_applycal = flagdata(vis=self.output_ms, mode='summary', name='applycal')
        status, report3 = th.check_val(flags_applycal['flagged'], self.exp_flags['flags_applycal'], exact=True, \
                               valname='flags_applycal')

        # Create the calibrator image and check the flux density and beam
        tclean(vis=self.output_ms,field='J0501-0159',imsize=256,cell='0.25arcsec',spw='16,18,20,22',
               imagename=self.image_prefix,niter=5000,nsigma=3.0)
        
        # Check calibrator peak flux density matches <1% of flux: 1.251 Jy 
        flux = imstat(imagename=self.image_prefix+'.image')['flux']
        status, report4 = th.check_val(flux, self.exp_flux, exact=False, epsilon=self.exp_tol, \
                               valname='flux')

        # Check that the beam maj/minor are within 1%; flux: 1.553 Jy; beam: 5.035" x 2.975" 78.692 deg
        beam = imhead(imagename=self.image_prefix+'.image')['restoringbeam']
        status, report5 = th.check_val(beam['major']['value'], self.exp_beam['major'], exact=False, 
                                       epsilon=self.exp_tol, valname='beam_major')

        status, report6 = th.check_val(beam['minor']['value'], self.exp_beam['minor'], exact=False, 
                                       epsilon=self.exp_tol, valname='beam_minor')
        
        status, report7 = th.check_val(beam['positionangle']['value'], self.exp_beam['positionangle'], 
                                       exact=False, epsilon=self.exp_tol, valname='positionangle') 

        # Concatenate the string returned by each report into a single string
        report = report0 + report1 + report2 + report3 + report4 + report5 + report6 + report7

        failed = self.filter_report(str(report))
        test_dict[test_name]['report'] = report
        add_to_dict(self, output = test_dict, dataset = self.output_ms)
        
        self.assertTrue(th.check_final(pstr = report), msg = failed)


##########   VLA test ############
class Test_srdp_vla(Test_srdp_base):
    """ VLA dataset """

    def setUp(self) -> None:
        # Reference criteria using casa-6.5.4-9-pipeline-2023.1.0.125
        # Expected flags and visibilities
        self.exp_total_vis = 1339392000.0
        self.exp_flags = {'flags_importasdm':0.0, 'flags_restored':629265468.0, 'flags_applycal':629265468.0,\
                          'flags_statwt':629290104.0}

        # Flux density in Jy/beam; major, minor in arcsec and position angle in deg
        self.exp_statwt = {'mean':115.30387726818151,"variance":7210.116933975743}
        self.exp_image = {'flux':0.7645,'major':21.292,'minor':9.361,"positionangle":1.405}
        self.exp_tol = 0.01 # within 1%

        # Input data
        self.asdm_name = '18A-426.sb35644955.eb35676220.58411.96917952546'
        self.restore_flags_tar = self.asdm_name + '.ms.flagversions.tgz'
        self.caltables = ['18A-426.sb35644955.eb35676220.58411.96917952546.ms.hifv_priorcals.s5_3.gc.tbl',
                          '18A-426.sb35644955.eb35676220.58411.96917952546.ms.hifv_priorcals.s5_4.opac.tbl',
                          '18A-426.sb35644955.eb35676220.58411.96917952546.ms.hifv_priorcals.s5_5.rq.tbl',
                          '18A-426.sb35644955.eb35676220.58411.96917952546.ms.finaldelay.k',
                          '18A-426.sb35644955.eb35676220.58411.96917952546.ms.finalBPcal.b',
                          '18A-426.sb35644955.eb35676220.58411.96917952546.ms.averagephasegain.g',
                          '18A-426.sb35644955.eb35676220.58411.96917952546.ms.finalampgaincal.g',
                          '18A-426.sb35644955.eb35676220.58411.96917952546.ms.finalphasegaincal.g']

        # Get input data
        os.symlink(data_path + self.asdm_name, self.asdm_name)
        os.symlink(data_path + self.restore_flags_tar, self.restore_flags_tar)

        # Copy cal tables locally
        for ctable in self.caltables:
            shutil.copytree(data_path + ctable, ctable)

        # Output data
        self.output_ms = self.asdm_name + '.ms'
        self.restore_flags = self.output_ms + '.flagversions'
        self.image_prefix = 'J1820-2528'


    def tearDown(self) -> None:
        os.system('rm -rf '+self.image_prefix + '*')
        shutil.rmtree(self.output_ms)

        os.unlink(self.asdm_name)
        os.unlink(self.restore_flags_tar)
        shutil.rmtree(self.restore_flags)
        for ctable in self.caltables:
            shutil.rmtree(ctable)

    @stats_dict(test_dict)  
    def test_vla_restore(self):
        """ SRDP VLA pipeline restore mode: check that restored flags do not change """
        
        test_name = self._testMethodName
        report = []
        # Import the MS and do not apply the online flags
        # TODO: check the outfile usage
        importasdm(asdm=self.asdm_name, vis=self.output_ms, createmms=False, ocorr_mode='co',
                   lazy=False, asis='Receiver CalAtmosphere', process_caldevice=True,
                   process_pointing=True, savecmds=False, overwrite=False, bdfflags=False, 
                   with_pointing_correction=True)

        # There should be zero flags in the created MS
        flags_importasdm = flagdata(vis=self.output_ms, mode='summary', name='after importasdm')

        # Check the total visibilities
        status, report0 = th.check_val(flags_importasdm['total'], self.exp_total_vis, exact=True, \
                               valname='total_visibilities')

        # Check the flags after importasdm
        status, report1 = th.check_val(flags_importasdm['flagged'], self.exp_flags['flags_importasdm'], exact=True, \
                               valname='flags_importasdm')

        # Run hanningsmooth and create a new MS. Rename new one to original MS
        hanningsmooth(vis=self.output_ms, outputvis='temphanning.ms', datacolumn='data')
        shutil.rmtree(self.output_ms, ignore_errors=False)
        shutil.move('temphanning.ms', self.output_ms)

        # Restore the Pipeline_Final flags to the MS and check flags
        os.system('tar xvf 18A-426.sb35644955.eb35676220.58411.96917952546.ms.flagversions.tgz')
        flagmanager(vis=self.output_ms, mode='restore', versionname='Pipeline_Final')
        flags_restored = flagdata(vis=self.output_ms, mode='summary', name='before-applycal')

        # Check the flags before applycal, which should contain only the restored flags from flagmanager
        status, report2 = th.check_val(flags_restored['flagged'], self.exp_flags['flags_restored'], exact=True, \
                               valname='flags_restored')

        # Run applycal
        applycal(vis=self.output_ms,antenna='*&*', gaintable=self.caltables, gainfield=['', '', '', '', '', '', '', ''],
                 interp=['linear','linear','linear','linear','linear,nearestflag','linear','linear','linear'],
                 spwmap=[[], [], [], [], [], [], [], []], calwt=[False,False,False,False,False,False,False,False],
                 parang=True, applymode='calflagstrict', flagbackup=False)

        # Check that the flags did not change
        flags_applycal = flagdata(vis=self.output_ms, mode='summary', name='after applycal')
        status, report3 = th.check_val(flags_applycal['flagged'], self.exp_flags['flags_applycal'], exact=True, \
                               valname='flags_applycal')

        # Get weigth statistics
        # TODO: get a tolerance for the statwt comparisons
        stats = statwt(vis=self.output_ms, minsamp=8, datacolumn='corrected')

        # These values are using the default th.check_val tolerance of 0.05
        status, report4 = th.check_val(stats['mean'], self.exp_statwt['mean'], exact=False, valname='statwt_mean')
        status, report5 = th.check_val(stats['variance'], self.exp_statwt['variance'], exact=False, valname='statwt_variance')

        # Check the statwt flags
        flags_statwt = flagdata(vis=self.output_ms, mode='summary', name='after statwt')
        status, report6 = th.check_val(flags_statwt['flagged'], self.exp_flags['flags_statwt'], exact=False, \
                               valname='flags_statwt')

        # Create calibrator image
        tclean(vis=self.output_ms, field='J1820-2528', imsize=256, cell='1.0arcsec', spw='0~15',
        imagename=self.image_prefix, niter=5000, nsigma=5.0)

        # Calibrator flux density and beam: 0.7645 Jy/beam, 21.292" x 9.361" 1.405 degrees
        flux = imstat(imagename=self.image_prefix+'.image')['flux']
        status, report7 = th.check_val(flux, self.exp_image['flux'], exact=False, epsilon=self.exp_tol, \
                               valname='flux')
        beam = imhead(imagename=self.image_prefix+'.image')['restoringbeam']
        status, report8 = th.check_val(beam['major']['value'], self.exp_image['major'], exact=False, epsilon=self.exp_tol, \
                               valname='beam_major')
        status, report9 = th.check_val(beam['minor']['value'], self.exp_image['minor'], exact=False, epsilon=self.exp_tol, \
                               valname='beam_minor')
        status, report10 = th.check_val(beam['positionangle']['value'], self.exp_image['positionangle'], exact=False, \
                                epsilon=self.exp_tol, valname='positionangle')
       
        # Concatenate the string returned by each report into a single string
        # Each report# is a tuple(Bool,string)
        report = report0 + report1 + report2 + report3 + report4 + report5 + report6 + \
                    report7 + report8 + report9 + report10
        
        failed = self.filter_report(str(report))
        test_dict[test_name]['report'] = report
        add_to_dict(self, output = test_dict, dataset = self.output_ms)
        
        self.assertTrue(th.check_final(pstr = report), msg = failed)

if __name__ == '__main__':
    unittest.main()
