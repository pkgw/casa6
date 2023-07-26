##########################################################################
# test_task_phaseshift.py
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
# [Add the link to the JIRA ticket here once it exists]
#
# Based on the requirements listed in plone found here:
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.manipulation.phaseshift.html
#
#
##########################################################################
import glob
import numpy as np
import os
import shutil
import unittest

from casatools import componentlist, measures

from casatasks import calmod

class calmod_test(unittest.TestCase):


    def setUp(self):
        self.cl = componentlist()


    def tearDown(self):
        self.cl.done()
        del self.cl

    
    def exception_verification(self, cm, expected_msg):
        exc = cm.exception
        pos = str(exc).find(expected_msg)
        self.assertNotEqual(
            pos, -1, msg=f'Unexpected exception was thrown: {exc}'
        )


    def test_inputs(self):
        '''Test inputs meet various constraints'''
        with self.assertRaises(ValueError) as cm: 
            calmod()
        self.exception_verification(cm, 'outfile must be specified')
        with self.assertRaises(ValueError) as cm: 
            calmod('my.cl')
        self.exception_verification(cm, 'Exactly one of source or direction must be specified')
        with self.assertRaises(ValueError) as cm: 
            calmod('my.cl', 'mysource', 'mydirection')
        self.exception_verification(cm, 'Both source and direction may not be simultaneously specified')
        with self.assertRaises(ValueError) as cm: 
            calmod('my.cl', 'mysource')
        self.exception_verification(cm, 'Unsupported calibrator mysource')
        with self.assertRaises(ValueError) as cm:
            calmod('my.cl', direction='mydirection')
        self.exception_verification(cm, 'Illegal direction specification mydirection')
        with self.assertRaises(ValueError) as cm:
            calmod('my.cl', direction='1 2 3')
        self.exception_verification(cm, 'Illegal direction specification 1 2 3')
        with self.assertRaises(ValueError) as cm: 
            calmod('my.cl', '3c48')
        self.exception_verification(cm, 'band must be specified')
        with self.assertRaises(ValueError) as cm: 
            calmod('my.cl', '3c48', band='m')
        self.exception_verification(cm, 'band m not supported')
        with self.assertRaises(ValueError) as cm: 
            calmod('my.cl', '3c48', band='q', obsdate=0)
        self.exception_verification(cm, 'obsdate must be >= 44239')
        with self.assertRaises(ValueError) as cm: 
            calmod('my.cl', '3c48', band='q', obsdate=50000, refdate=1)
        self.exception_verification(cm, 'refdate must be <= 0 or >= 44239')
        with self.assertRaises(ValueError) as cm: 
            calmod('my.cl', '3c48', band='q', obsdate=50000, refdate=0)
        self.exception_verification(cm, 'hosts must be specified')
        hosts = ['https://my.out.edu', 'zz']
        with self.assertRaises(ValueError) as cm: 
            calmod('my.cl', '3c48', band='q', obsdate=50000, refdate=0, hosts=hosts)
        self.exception_verification(cm, 'zz is not a valid host expressed as a URL')
        

if __name__ == '__main__':
     unittest.main()
