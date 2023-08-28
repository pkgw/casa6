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
import subprocess
import sys
import threading
from time import sleep
import unittest
from urllib import request
from urllib.error import URLError


from casatasks import casalog

from casatools import componentlist, measures

from casatasks import calmod

import casatestutils

class calmod_test(unittest.TestCase):


    hostname = 'http://127.0.0.1:8080'


    @classmethod
    def capture_output(cls, subprocess):
        stdout, stderr = cls.web_server.communicate()
        casalog.post(f'stdout: {stdout.decode("utf-8")}', 'INFO')
        casalog.post(f'stderr: {stderr.decode("utf-8")}', 'WARN')



    @classmethod
    def setUpClass(cls):

        server = os.sep.join([casatestutils.__path__[0],
            'calmod_helpers', 'vla_mock_server.py'])
        casalog.post(f'server is {server}', 'INFO')
        cls.web_server = subprocess.Popen(
            [sys.executable, server], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        url = 'http://127.0.0.1:8080'
        req = request.Request(url)
        i = 0
        started = False
        while not started:
            try:
                with request.urlopen(req) as response:
                    pass
            except (ConnectionRefusedError, URLError) as e:
                started = str(e) == 'HTTP Error 400: BAD REQUEST'
                if i == 19:
                    raise RuntimeError('Unable to start web server within 10 seconds')
            i += 1
            sleep(0.5)
            if started:
                break
        casalog.post(
            f'Web server successfully started between {0.5*(i-1)} and {0.5*i} seconds',
            'INFO'
        )


    def setUp(self):
        self.cl = componentlist()
        self.clname = 'my.cl'


    def tearDown(self):
        self.cl.done()
        del self.cl
        if os.path.exists(self.clname):
            shutil.rmtree(self.clname)

    
    @classmethod
    def tearDownClass(cls):
        output_thread = threading.Thread(
            target=cls.capture_output, args=(cls.web_server,)
        )
        output_thread.start()
        sleep(2)
        cls.web_server.terminate()
        output_thread.join()


    def exception_verification(self, cm, expected_msg):
        exc = cm.exception
        pos = str(exc).find(expected_msg)
        self.assertNotEqual(
            pos, -1, msg=f'Unexpected exception was thrown: {exc}'
        )


    def test_inputs(self):
        """Test inputs meet various constraints"""
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
            calmod('my.cl', '3c48', band='q', obsdate=1)
        self.exception_verification(cm, 'obsdate must be <= 0 or >= 44239')
        with self.assertRaises(ValueError) as cm: 
            calmod('my.cl', '3c48', band='q', obsdate=50000, refdate=1)
        self.exception_verification(cm, 'refdate must be <= 0 or >= 44239')
        with self.assertRaises(ValueError) as cm: 
            calmod('my.cl', '3c48', band='q', obsdate=50000, refdate=0, hosts=[])
        self.exception_verification(cm, 'hosts must be specified')
        hosts = ['zz']
        with self.assertRaises(ValueError) as cm: 
            calmod('my.cl', '3c48', band='q', obsdate=50000, refdate=0, hosts=hosts)
        self.exception_verification(cm, 'zz is not a valid host expressed as a URL')
        hosts = ['http://my.bogus.com:8080']
        with self.assertRaises(RuntimeError) as cm: 
            calmod('my.cl', '3c48', band='q', obsdate=50000, refdate=0, hosts=hosts)
        self.exception_verification(cm, 'All URLs failed to return a component list')


    def test_component_list_writing(self):
        """Test successful writing of a component list"""
        hosts = [self.hostname]
        calmod(self.clname, '3C48', band='Q',obsdate=50000, hosts=hosts)
        self.cl.open(self.clname)
        self.assertEqual(self.cl.length(), 385, 'Incorrect number of components')
        ws = self.cl.getkeyword('web_service')
        self.assertEqual(ws['band'], 'Q', 'Incorrect band in web_service metadata')
        self.assertEqual(ws['source'], '3C48', 'Incorrect source in web_service metadata')
        
    
    def test_direction(self):
        """Test direction input"""
        hosts = [self.hostname]
        direction = 'J2000 01:37:41.1 33.09.32'
        calmod(self.clname, direction=direction, band='Q',obsdate=50000, hosts=hosts)
        self.cl.open(self.clname)
        self.assertEqual(self.cl.length(), 385, 'Incorrect number of components')
        ws = self.cl.getkeyword('web_service')
        self.assertEqual(ws['band'], 'Q', 'Incorrect band in web_service metadata')


if __name__ == '__main__':
     unittest.main()
