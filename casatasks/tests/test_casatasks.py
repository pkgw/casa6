#########################################################################
# test_casatasks.py
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
#
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.html
#
##########################################################################
import os
import unittest
import numpy as np
from casatools import table
from casatasks import importasdm, importfits, listobs, flagdata, gencal, setjy, fluxscale
from casatasks import gaincal
from casatestutils import sparse_check


class BaseClass(unittest.TestCase):
    def getdata(testfiles=None):
        # Download data for the tests
        sparse_check.download_data(testfiles)

class CasaTasksTests(BaseClass):
    """ Unit tests for CASA tasks
        The test script will download the test data automatically to
        the local directory and will remove them afterwards
        """
    @classmethod
    def setUpClass(cls) -> None:
        cls.asdm = 'AutocorrASDM'
        cls.fitsimage = 'two_gaussian_model.fits'
        cls.listobs_ms = 'uid___X02_X3d737_X1_01_small.ms'
        cls.dummy = 'pm_ngc5921.ms'
        cls.flagdata_ms = 'ngc5921.ms'
        cls.gencal_ms = 'tdem0003gencal.ms'
        cls.setjy_ms = 'ngc5921.ms'
        cls.fluxscale_ms = 'CalMSwithModel.ms'
        cls.fluxscale_gtable = 'ModelGcal.G0'
        cls.gaincal_ms = 'gaincaltest2.ms'
        cls.input_files = [cls.asdm, cls.fitsimage, cls.listobs_ms, cls.dummy, cls.flagdata_ms,
                           cls.gencal_ms, cls.fluxscale_ms, cls.fluxscale_gtable, cls.gaincal_ms]
        # Fetch input data
        cls.getdata(testfiles=cls.input_files)

        # Output data
        cls.asdm_ms = cls.asdm + '.ms'
        cls.onlineflags = 'onlineflags.txt'
        cls.casaimage = cls.fitsimage + '.image'
        cls.gentable = 'gencal_antpos.cal'
        cls.fluxscale_out = 'fluxout.cal'
        cls.gaincal_out = 'gaincaltable.cal'
        cls.output_files = [cls.asdm_ms, cls.onlineflags, cls.casaimage, cls.gentable, cls.fluxscale_out,
                            cls.gaincal_out]
    @classmethod
    def tearDownClass(cls) -> None:
        # Remove input files
        for inpfile in cls.input_files:
            os.system('rm -rf ' + inpfile)

        # Remove output files
        for outfile in cls.output_files:
            os.system('rm -rf ' + outfile)

    def setUp(self):
        print("%s: %s" % (self._testMethodName, self._testMethodDoc))

    def test_importasdm_autocorrelations(self):
        """Test importasdm on autocorrelation ASDM with scan selection and saving online flags"""
        # Use default name for output MS
        self.assertEqual(
            importasdm(asdm=self.asdm, scans='3', savecmds=True, outfile=self.onlineflags, flagbackup=False), None)
        self.assertTrue(os.path.exists(self.onlineflags))
        with open(self.onlineflags, 'r') as ff:
            cmds = ff.readlines()
            # auto-correlation should have been written to online flags
            self.assertTrue(cmds[0].__contains__('&&*'))

    def test_importfits_beam(self):
        """Test importfits on data with synthesized beam"""
        importfits(fitsimage=self.fitsimage, imagename=self.casaimage, beam=['0.35arcsec', '0.24arcsec', '25deg'])
        self.assertTrue(os.path.exists(self.casaimage), "Output image does not exist")

    def test_listobs_field_selection(self):
        """Test listobs on MS with field selection"""
        res = listobs(vis=self.listobs_ms, field='1')
        self.assertEqual(res['nfields'], 1)

    def test_flagdata_manual(self):
        """Test flagdata with manual mode"""
        flagdata(vis=self.flagdata_ms, correlation='LL', savepars=False, flagbackup=False)
        flagdata(vis=self.flagdata_ms, spw='0:17~19', savepars=False, flagbackup=False)
        flagdata(vis=self.flagdata_ms, antenna='VA05&&VA09', savepars=False, flagbackup=False)
        flagdata(vis=self.flagdata_ms, field='1', savepars=False, flagbackup=False)
        s = flagdata(vis=self.flagdata_ms, mode='summary', minrel=0.9, spwchan=True, basecnt=True)
        assert 'VA05&&VA09' in s['baseline']
        assert set(s['spw:channel']) == set(['0:17', '0:18', '0:19'])
        assert list(s['correlation'].keys()) == ['LL']  # LL
        assert list(s['field'].keys()) == ['1445+09900002_0']
        assert set(s['scan']) == set(['2', '4', '5', '7'])  # field 1
        s = flagdata(vis=self.flagdata_ms, mode='summary', maxrel=0.8)
        assert set(s['field']) == set(['1331+30500002_0', 'N5921_2'])
        s = flagdata(vis=self.flagdata_ms, mode='summary', minabs=400000)
        assert set(s['scan']) == set(['3', '6'])
        s = flagdata(vis=self.flagdata_ms, mode='summary', minabs=400000, maxabs=450000)
        assert list(s['scan'].keys()) == ['3']

    def test_flagdata_quack(self):
        """Test flagdata with quack mode"""
        flagdata(vis=self.flagdata_ms, mode='unflag', flagbackup=False)
        flagdata(vis=self.flagdata_ms, mode='quack', quackmode='beg', quackinterval=1, flagbackup=False)
        self.assertEqual(flagdata(vis=self.flagdata_ms, mode='summary')['flagged'], 329994)

    def test_gencal_antpos_manual(self):
        """Test gencal with manual antenna position correction"""
        gencal(vis=self.gencal_ms, caltable=self.gentable, caltype='antpos',
               antenna='ea12,ea22', parameter=[-0.0072, 0.0045, -0.0017, -0.0220, 0.0040, -0.0190])
        self.assertTrue(os.path.exists(self.gentable))

    def test_setjy_flux_standard(self):
        """Test setjy with input dictionary for fluxscale standard"""
        fluxscaledict = {
            '1':
                {'0':
                     {'fluxd': np.array([2.48362403, 0., 0., 0.]),
                         'fluxdErr': np.array([0.00215907, 0., 0., 0.]),
                         'numSol': np.array([54., 0., 0., 0.])
                      },
                   'fieldName': '1445+09900002_0',
                   'fitFluxd': 0.0,
                   'fitFluxdErr': 0.0,
                   'fitRefFreq': 0.0,
                   'spidx': np.array([0., 0., 0.]),
                   'spidxerr': np.array([0., 0., 0.])
                 },
             'freq': np.array([1.41266507e+09]),
             'spwID': np.array([0], dtype=np.int32),
             'spwName': np.array(['none'], dtype='|S5')
        }

        retdict = setjy(vis=self.setjy_ms, standard='fluxscale',fluxdict=fluxscaledict,usescratch=False)
        self.assertEqual(retdict['1']['0']['fluxd'][0],2.48362403,0.0001)

    def test_fluxscale_refspwmap(self):
        """Test fluxscale with a refspwmap parameter"""
        fluxscale(vis=self.fluxscale_ms, caltable=self.fluxscale_gtable, fluxtable=self.fluxscale_out,
                  reference=['0'], refspwmap=[1, 1, 1, 1])
        self.assertTrue(os.path.exists(self.fluxscale_out))

    def test_gaincal_gaintype_g(self):
        """Test gaincal using gaintype G and flagged antennas"""
        gaincal(vis=self.gaincal_ms, caltable=self.gaincal_out, refant='0', field='0', solint='inf',
                combine='scan', antenna='0~5&', smodel=[1,0,0,0], gaintype='G')
        self.assertTrue(os.path.exists(self.gaincal_out))

if __name__ == '__main__':
    unittest.main()
