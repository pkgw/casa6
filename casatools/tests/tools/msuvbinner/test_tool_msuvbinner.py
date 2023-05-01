#########################################################################
# test_tool_msuvbinner.py
# Copyright (C) 2022
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatools.msuvbinner.html
#
# Testing of methods selectdata, setoutputms, filloutputms
#
##########################################################################

import os
import sys
import unittest
import numpy as np
import shutil
from casatools import ctsys, msuvbinner, table
ctsys_resolve = ctsys.resolve
from casatestutils import testhelper as th
datapath = ctsys_resolve('unittest/msuvbinner/')
_tb=table()
class msuvbinnerTest(unittest.TestCase):
    ms1='refim_point.ms'
    outms='tst_uvgrid.ms'
    def setUp(self):
         if(not os.path.exists(self.ms1)):
             os.system('cp -RH '+datapath+'/'+self.ms1+' .')
    def tearDown(self):
        shutil.rmtree(self.ms1)
        shutil.rmtree(self.outms, ignore_errors=True)
    def test_fill_1ms(self):
        '''test uvgridding 1 ms onto an output grid '''
        msbinner=msuvbinner.msuvbinner(phasecenter='J2000 19:59:28.5 40d44m01.5', nx=100, ny=100, ncorr=1, nchan=20, cellx='8arcsec', celly='8arcsec', fstart='975MHz', fstep='50MHz', memfrac=0.5, wproject=False, doflag=False)
        msbinner.selectdata(msname=self.ms1, spw='', field='0', taql='')
        msbinner.setoutputms(self.outms)
        msbinner.filloutputms()
        del msbinner
        _tb.open(self.outms)
        arr=_tb.getcol('DATA')
        _tb.done()
        self.assertAlmostEqual(np.max(arr[0,0,:]), 1.5, 4)
        self.assertAlmostEqual(np.max(arr[0,19,:]), 0.769271, 5)


if __name__ == '__main__':
    unittest.main()
 
