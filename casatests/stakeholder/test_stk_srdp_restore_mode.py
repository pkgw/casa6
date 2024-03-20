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
# # https://open-jira.nrao.edu/browse/CAS-14319
#
##########################################################################

'''
These tests will run importadsm to fill a CASA MS and restore the flag backup saved
in a previoun run of the same project. The tests will then check the flagging stats
between running several calibration tasks to verify that the flags have not
changed compared to the expected values from the repvious run. The last step will
create a calibrator image and check the flux density and beam values.
'''

import os
import shutil
import re
import glob
import unittest
from casatools import ctsys
from casatasks import importasdm, flagdata, flagmanager, applycal, tclean, imstat, imhead, hanningsmooth, statwt

data_path = ctsys.resolve('stakeholder/srdp/')

def _extract_caltables(filename):
    # This pattern matches the 'caltable' key and extracts its value
    pattern = re.compile(r"caltable='(.*?)'")
    caltables = []

    with open(filename, 'r') as file:
        for line in file:
            match = pattern.search(line)
            if match:
                # Extract the value of caltable and add it to the list
                caltables.append(match.group(1))
        unique_tables = set(caltables)

    return unique_tables


class Test_srdp_alma_12m(unittest.TestCase):
    """ ALMA 12m dataset """

    def setUp(self) -> None:
        # Reference criteria using casa-6.5.4-9-pipeline-2023.1.0.125
        # flux density in Jy; major, minor in arcsec and position angle in deg
        self.ref_flux = 1.251
        self.ref_beam = {'major':0.543,'minor':0.345,"positionangle":-19.153}
        self.tolerance = 0.01 # within 1%

        # Input data
        self.asdm_name = 'uid___A002_Xd3607d_X2b0d'
        self.restore_flags_tar = 'uid___A002_Xd3607d_X2b0d.ms.flagversions.tgz'
        # Reference callibrary has the names of the cal tables to apply
        self.reference_callibrary = 'uid___A002_Xd3607d_X2b0d.ms.s1.3.callibrary'

        self.output_ms = self.asdm_name + '.ms'
        self.restore_flags = 'uid___A002_Xd3607d_X2b0d.ms.flagversions'
        self.image_prefix = 'J0336+3218'

        # Get input data
        os.symlink(data_path + self.asdm_name, self.asdm_name)
        os.symlink(data_path + self.restore_flags_tar, self.restore_flags_tar)
        shutil.copy(data_path + self.reference_callibrary, self.reference_callibrary)
        self.caltables = _extract_caltables(self.reference_callibrary)

        # Copy cal tables locally
        for ctable in self.caltables:
            shutil.copytree(data_path + ctable, ctable)

    def tearDown(self) -> None:
        os.system('rm -rf '+self.image_prefix + '*')
        shutil.rmtree(self.output_ms)

        os.unlink(self.asdm_name)
        os.unlink(self.restore_flags_tar)
        shutil.rmtree(self.restore_flags)
        os.remove(self.reference_callibrary)
        for ctable in self.caltables:
            shutil.rmtree(ctable)

    def test_alma_12m_restore(self):
        """SRDP stakeholder alma 12m: Check that restored flags do not change in new MS"""

        # Import the MS and do not apply the online flags
        # The BDF flags are saved with the flagbackup in the .flagversions file with name Original
        importasdm(asdm=self.asdm_name, vis=self.output_ms, createmms=False, ocorr_mode='ca', lazy=False,
                   asis='SBSummary ExecBlock Antenna Station Receiver Source CalAtmosphere CalWVR CalPointing',
                   process_caldevice=False, process_flags=True, applyflags=False, savecmds=False,
                   overwrite=False, bdfflags=True, with_pointing_correction=False)

        # Get the flags summary and compare with the expected values
        flags_before = flagdata(vis=self.output_ms, mode='summary', name='before-flagmanager')
        self.assertEqual(flags_before['flagged'], 3900108.0, 'Flagged visibilities do not match')
        self.assertEqual(flags_before['total'], 1092709620.0, 'Total visibilities do not match')

        # Restore the reference flags from the pipeline to the MS
        os.system('tar xf '+ self.restore_flags_tar)
        flagmanager(vis=self.output_ms, mode='restore', versionname='Pipeline_Final')

        # Get the flags summary and compare with expected values
        # At this point, the MS has the sum of flags from Original and Pipeline_Final
        flags_after = flagdata(vis='uid___A002_Xd3607d_X2b0d.ms', mode='summary', name='before-applycal')
        self.assertEqual(flags_after['flagged'], 159843772.0, 'Flagged visibilities do not match')
        self.assertEqual(flags_after['total'], 1092709620.0, 'Total visibilities do not match')

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
        self.assertEqual(flags_final['flagged'], 159843772.0, 'Flagged visibilities do not match')
        self.assertEqual(flags_final['total'], 1092709620.0, 'Total visibilities do not match')

        # Create the calibrator image and check the flux density and beam
        tclean(vis=self.output_ms,field='J0336+3218',imsize=256,cell='0.05arcsec',spw='19,21,23',
               imagename=self.image_prefix,niter=5000,nsigma=3.0)

        # Check calibrator peak flux density matches <1%, beam maj/minor within 1%; flux: 1.251 Jy; beam: 0.543" x 0.345" -19.153 deg
        flux = imstat(imagename=self.image_prefix+'.image')['flux']
        beam = imhead(imagename=self.image_prefix+'.image')['restoringbeam']
        self.assertAlmostEqual(flux, self.ref_flux, delta=self.tolerance,
                               msg="Flux density does not match")
        self.assertAlmostEqual(beam['major']['value'],self.ref_beam['major'], delta=self.tolerance,
                               msg="Beam major does not match")
        self.assertAlmostEqual(beam['minor']['value'],self.ref_beam['minor'], delta=self.tolerance, 
                               msg="Beam minor does not match")
        self.assertAlmostEqual(beam['positionangle']['value'],self.ref_beam['positionangle'], delta=self.tolerance, 
                               msg="PA does not match")

class Test_srdp_vla(unittest.TestCase):
    """ VLA dataset """

    def setUp(self) -> None:
        # Reference criteria using casa-6.5.4-9-pipeline-2023.1.0.125
        # flux density in Jy/beam; major, minor in arcsec and position angle in deg
        self.ref_image = {'flux':0.7645,'major':21.292,'minor':9.361,"positionangle":1.405}
        self.ref_statwt = {'mean':115.30387726818151,"variance":7210.116933975743}
        self.tolerance = 0.01 # within 1%

        # I/O data
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

        self.output_ms = self.asdm_name + '.ms'
        self.restore_flags = self.output_ms + '.flagversions'
        self.image_prefix = 'J1820-2528'

        # Get input data
        os.symlink(data_path + self.asdm_name, self.asdm_name)
        os.symlink(data_path + self.restore_flags_tar, self.restore_flags_tar)

        # Copy cal tables locally
        for ctable in self.caltables:
            shutil.copytree(data_path + ctable, ctable)

    def tearDown(self) -> None:
        #os.system('rm -rf '+self.image_prefix + '*')
        shutil.rmtree(self.output_ms)

        os.unlink(self.asdm_name)
        os.unlink(self.restore_flags_tar)
        shutil.rmtree(self.restore_flags)
        for ctable in self.caltables:
            shutil.rmtree(ctable)

    def test_vla_restore(self):
        """SRDP stakeholder vla: Check that restored flags do not change"""

        # Import the MS and do not apply the online flags
        # TODO: check the outfile usage
        importasdm(asdm=self.asdm_name, vis=self.output_ms, createmms=False, ocorr_mode='co',
                   lazy=False, asis='Receiver CalAtmosphere', process_caldevice=True,
                   process_pointing=True, savecmds=False, overwrite=False, bdfflags=False, 
                   with_pointing_correction=True)

        # There should be zero flags in the created MS
        flags_importasdm = flagdata(vis=self.output_ms, mode='summary', name='after importasdm')
        self.assertEqual(flags_importasdm['flagged'], 0, 'Flagged visibilities shuold be zero')
        self.assertEqual(flags_importasdm['total'],1339392000.0, 'Total visibilities do not match')

        # Run hanningsmooth and create a new MS. Rename new one to original MS
        hanningsmooth(vis=self.output_ms, outputvis='temphanning.ms', datacolumn='data')
        shutil.rmtree(self.output_ms, ignore_errors=False)
        shutil.move('temphanning.ms', self.output_ms)

        # Restore the Pipeline_Final flags to the MS and check flags
        os.system('tar xvf 18A-426.sb35644955.eb35676220.58411.96917952546.ms.flagversions.tgz')
        flagmanager(vis=self.output_ms, mode='restore', versionname='Pipeline_Final')
        flags_before_applycal = flagdata(vis=self.output_ms, mode='summary', name='before-applycal')
        self.assertEqual(flags_before_applycal['flagged'], 629265468.0, 'Flagged visibilities do not match')
        self.assertEqual(flags_before_applycal['total'],1339392000.0, 'Total visibilities do not match')

        # Run applycal
        applycal(vis=self.output_ms,antenna='*&*', gaintable=self.caltables, gainfield=['', '', '', '', '', '', '', ''],
                 interp=['linear','linear','linear','linear','linear,nearestflag','linear','linear','linear'],
                 spwmap=[[], [], [], [], [], [], [], []], calwt=[False,False,False,False,False,False,False,False],
                 parang=True, applymode='calflagstrict', flagbackup=False)

        # Check that the flags did not change
        flags_after_applycal = flagdata(vis=self.output_ms, mode='summary', name='after applycal')
        self.assertEqual(flags_after_applycal['flagged'], 629265468.0, 'Flagged visibilities do not match')
        self.assertEqual(flags_after_applycal['total'],1339392000.0, 'Total visibilities do not match')

        # Get weigth statistics
        # TODO: get a tolerance for these comparisons
        stats = statwt(vis=self.output_ms, minsamp=8, datacolumn='corrected')
        self.assertAlmostEqual(stats['mean'], self.ref_statwt['mean'], delta=self.tolerance,
                         msg="Statwt mean does not match")
        self.assertAlmostEqual(stats['variance'], self.ref_statwt['variance'], delta=self.tolerance,
                         msg="Statwt variance does not match")

        # Check the flags
        flags_statwt = flagdata(vis=self.output_ms, mode='summary', name='after statwt')
        self.assertAlmostEqual(flags_statwt['flagged'], 629290104.0, 'Flagged visibilities do not match')
        self.assertAlmostEqual(flags_statwt['total'],1339392000.0, 'Total visibilities do not match')

        # Create calibrator image
        tclean(vis=self.output_ms, field='J1820-2528', imsize=256, cell='1.0arcsec', spw='0~15',
        imagename=self.image_prefix, niter=5000, nsigma=5.0)

        # Calibrator flux density and beam: 0.7645 Jy/beam, 21.292" x 9.361" 1.405 degrees
        flux = imstat(imagename=self.image_prefix+'.image')['flux']
        beam = imhead(imagename=self.image_prefix+'.image')['restoringbeam']
        self.assertAlmostEqual(flux, self.ref_image['flux'], delta=self.tolerance,
                               msg="Flux density does not match")
        self.assertAlmostEqual(beam['major']['value'],self.ref_image['major'], delta=self.tolerance,
                               msg="Beam major does not match")
        self.assertAlmostEqual(beam['minor']['value'],self.ref_image['minor'], delta=self.tolerance,
                               msg="Beam minor does not match")
        self.assertAlmostEqual(beam['positionangle']['value'],self.ref_image['positionangle'], delta=self.tolerance,
                               msg="PA does not match")
        


if __name__ == '__main__':
    unittest.main()
