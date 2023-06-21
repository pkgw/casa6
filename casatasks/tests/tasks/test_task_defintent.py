#########################################################################
# test_task_defintent.py
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.calibration.defintent.html
#
#
##########################################################################
import os
import sys
import shutil
import unittest
import itertools
import numpy as np

from casatools import ms, ctsys, table
from casatasks import defintent

ctsys_resolve = ctsys.resolve

tb = table()


datapath = ctsys.resolve('unittest/defintent/')

class Defintent_tests(unittest.TestCase):
    msfile = 'gaincaltest2.ms'

    def setUp(self):
        shutil.copytree(os.path.join(datapath,self.msfile), self.msfile)

    def tearDown(self):
        shutil.rmtree(self.msfile)
        
    def test_mode_append(self):
        """ A new intent should be added that appends to the previous intent """
        defintent(vis=self.msfile, intent='testintent', mode='append')
        
        # Get the intents and states
        tb.open(self.msfile + '/STATE')
        intents = tb.getcol('OBS_MODE')
        tb.close()
        tb.open(self.msfile)
        states = tb.getcol('STATE_ID')
        tb.close()
        
        # Check that the new intent was added to intents and that it hasn't overwritten
        self.assertTrue('testintent' in intents[-1])
        self.assertTrue('testintent' != intents[-1])
        # Check that the states were written to the new value
        self.assertTrue(np.all(states==1))
        
    def test_mode_set(self):
        """ A new intent is added and replaces the old intent """
        defintent(vis=self.msfile, intent='testintent', mode='set')
        
        # Get the intents and states
        tb.open(self.msfile + '/STATE')
        intents = tb.getcol('OBS_MODE')
        tb.close()
        tb.open(self.msfile)
        states = tb.getcol('STATE_ID')
        tb.close()
        
        # Check that the new intent was added to intents and that it has been
        self.assertTrue('testintent' == intents[-1])
        # Check that the states were written to the new value
        self.assertTrue(np.all(states == 1))
        
    def test_select_scan(self):
        """ A selection of intents to replace can be made with scan """
        # run defintent only selecting a subset of scans
        defintent(vis=self.msfile, intent='testintent', mode='set', scan='0,1,2,3')
        
        # Get the scans and new states
        tb.open(self.msfile)
        scanNumber = tb.getcol('SCAN_NUMBER')
        stateId = tb.getcol('STATE_ID')
        tb.close()
        
        # iterate through the rows and check the selected scans
        scans = [0,1,2,3]
        statesmatch = True
        for i in range(len(scanNumber)):
            if scanNumber[i] in scans and stateId[i] != 1:
                statesmatch = False
        # if the states weren't replaced for the scans then fail
        self.assertTrue(statesmatch)
        
    def test_select_field(self):
        """ A selection of intents to replace can be made with field """
        # run defintent only selecting a subset of fields
        defintent(vis=self.msfile, intent='testintent', mode='set', field='0')
        
        # Get the scans and new states
        tb.open(self.msfile)
        fieldNumber = tb.getcol('FIELD_ID')
        stateId = tb.getcol('STATE_ID')
        tb.close()
        
        # iterate through the rows and check the selected fields
        statesmatch = True
        for i in range(len(fieldNumber)):
            if fieldNumber[i] == 0 and stateId[i] != 1:
                statesmatch = False
        # if the states weren't replaced for the fields then fail
        self.assertTrue(statesmatch)
        
    def test_select_obsid(self):
        """ A selection of intents to replace can be made with obsid """
        # TODO: Look for ms with more interesting obs ids
        # run defintent only selecting a subset of obsids
        defintent(vis=self.msfile, intent='testintent', mode='set', obsid='0')
        
        # Get the scans and new states
        tb.open(self.msfile)
        obsNumber = tb.getcol('OBSERVATION_ID')
        stateId = tb.getcol('STATE_ID')
        tb.close()
        
        # iterate through the rows and check the selected obsids
        statesmatch = True
        for i in range(len(obsNumber)):
            if obsNumber[i] == 0 and stateId[i] != 1:
                statesmatch = False
        # if the states weren't replaced for the obsids then fail
        self.assertTrue(statesmatch)
    
if __name__ == '__main__':
    unittest.main()
