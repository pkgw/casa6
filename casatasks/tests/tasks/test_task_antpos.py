##########################################################################
# test_task_antpos.py
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
"""
import glob
import http.server
import numpy as np
import os
import re
import shutil
import sys
import threading
"""
import unittest
"""
from urllib import request
from urllib.error import URLError
from urllib.parse import urlparse, parse_qs

from casatasks import casalog

from casatools import componentlist, measures
"""

from casatasks import antpos

# import casatestutils

"""
# NOTE be certain to specify the top-level casatestutils directory
# in your PYTHONPATH so you load the casatestutils directory which
# is a subdir of that


class MockHTTPRequestHandler(http.server.BaseHTTPRequestHandler):
    """"HTTPServer mock request handler""""

    def do_GET(self):
        casalog.post('server path ' + self.path, 'WARN')
        parms = parse_qs(urlparse(self.path).query)
        casalog.post(f'server parms {parms}', 'INFO')
        good_sources = ("3C48", "3C286", "3C138", "3C147")
        if 'source' in parms and parms['source'][0].upper() not in good_sources:
            explain = f'source must be one of {good_sources}'
            self.send_error(400, message='Invalid input', explain=explain)
            self.end_headers()
            return
        """"Handle GET requests""""
        self.send_response(200)
        self.send_header("Content-Type", "application/json")
        self.end_headers()
        myfile = os.sep.join([
            casatestutils.__path__[0], 'calmod_helpers', 'query1.json'
        ])
        with open(myfile, 'r') as f:
            file_contents = f.read()
        self.wfile.write(str.encode(file_contents))

    def log_request(self, code=None, size=None):
        """"Don't log anything""""

"""

class antpos_test(unittest.TestCase):
    """

    hostname = 'http://127.0.0.1:8080'
    """

    def setUp(self):
        pass

    def tearDown(self):
        pass

    def exception_verification(self, cm, expected_msg):
        exc = cm.exception
        pos = str(exc).find(expected_msg)
        self.assertNotEqual(
            pos, -1, msg=f'Unexpected exception was thrown: {exc}'
        )

    """
    def query_server(self, method):
        server = http.server.ThreadingHTTPServer(
            ('127.0.0.1', 8080), MockHTTPRequestHandler
        )
        with server:
            server_thread = threading.Thread(target=server.serve_forever)
            server_thread.daemon = True
            server_thread.start()
            try:
                method()
            finally:
                server.shutdown()
    """

    def test_inputs(self):
        """Test inputs meet various constraints"""
        with self.assertRaises(ValueError) as cm: 
            antpos()
        self.exception_verification(cm, "Parameter outfile must be specified")
        with self.assertRaises(ValueError) as cm: 
            antpos(hosts=["good.example.com"])
        self.exception_verification(cm, "Parameter outfile must be specified")
        with self.assertRaises(ValueError) as cm: 
            antpos(outfile="myants.json", hosts=[])
        self.exception_verification(cm, "Parameter hosts must be specified")
        with self.assertRaises(ValueError) as cm: 
            antpos(outfile="myants.json", tw="1,2,3", hosts=["good.example.com"])
        self.exception_verification(
            cm, "Parameter tw should contain exactly one comma that separates two times"
        )
        with self.assertRaises(ValueError) as cm: 
            antpos(outfile="myants.json", tw="1,3", hosts=["good.example.com"])
        self.exception_verification(
            cm,
            "Begin time 1 does not appear to have a valid format. The correct format "
            "is of the form YYYY-MM-DDThh:mm:ss."
        )
        with self.assertRaises(ValueError) as cm: 
            antpos(
                outfile="myants.json", tw="2023-04-15T17:15:22,3",
                hosts=["good.example.com"]
            )
        self.exception_verification(
            cm,
            "End time 3 does not appear to have a valid format. The correct format "
            "is of the form YYYY-MM-DDThh:mm:ss."
        )
        with self.assertRaises(ValueError) as cm: 
            antpos(
                outfile="myants.json", tw="2023-04-15T17:15:22,2023-03-02T14:41:00",
                hosts=["good.example.com"]
            )
        self.exception_verification(
            cm,
            "Parameter tw, start time (2023-04-15T17:15:22) must be less than end "
            "time (2023-03-02T14:41:00)."
        )
        with self.assertRaises(ValueError) as cm: 
            antpos(
                outfile="myants.json", tw="2023-04-15T17:15:22,2023-04-15T17:15:22",
                hosts=["good.example.com"]
            )
        self.exception_verification(
            cm,
            "Parameter tw, start time (2023-04-15T17:15:22) must be less than end "
            "time (2023-04-15T17:15:22)."
        )
        with self.assertRaises(ValueError) as cm: 
            antpos(
                outfile="myants.json", snr=-1, hosts=["good.example.com"]
            )
        self.exception_verification(
            cm, "Parameter snr (-1.0) must be non-negative."
        )
        with self.assertRaises(ValueError) as cm: 
            antpos(
                outfile="myants.json", search="sr", hosts=["good.example.com"]
            )
        self.exception_verification(
            cm,
            "Parameter search (=sr) must have a value of either 'both_latest' "
            "or 'both_closest'."
        )
        with self.assertRaises(ValueError) as cm: 
            antpos(
                outfile="myants.json", hosts=["bogus./12?.example.com"]
            )
        self.exception_verification(
            cm,
            "Parameter hosts: bogus./12?.example.com is not a valid host expressed as a URL."
        )
        with self.assertRaises(RuntimeError) as cm: 
            antpos(
                outfile="myants.json", hosts=["http://www.bogus.edu"]
            )
        self.exception_verification(
            cm,
            "All URLs failed to return an antenna position list."
        )


    """
    def test_component_list_writing(self):
        """"Test successful writing of a component list""""
        hosts = [self.hostname]
        self.query_server(
            lambda: calmod(
                self.clname, '3C48', band='Q',obsdate=50000, hosts=hosts
            )
        )
        self.cl.open(self.clname)
        self.assertEqual(self.cl.length(), 385, 'Incorrect number of components')
        ws = self.cl.getkeyword('web_service')
        self.assertEqual(ws['band'], 'Q', 'Incorrect band in web_service metadata')
        self.assertEqual(ws['source'], '3C48', 'Incorrect source in web_service metadata')


    def test_bad_source_name(self):
        hosts = [self.hostname]
        with self.assertRaises(RuntimeError) as cm: 
            self.query_server(
                lambda: calmod(
                    'my.cl', 'mysource', band='L', hosts=[self.hostname],
                    obsdate=50000
                )
            )
        self.exception_verification(cm, 'All URLs failed to return a component list')
        found = False
        pattern = "source must be one of \('3C48', '3C286', '3C138', '3C147'\)"
        with open(casalog.logfile()) as logfile:
            for line in logfile:
                if re.search(pattern, line):
                    found = True
                    break
        self.assertTrue(found)


    def test_obsdate_as_string(self):
        hosts = [self.hostname]
        self.query_server(
            lambda: calmod(
                self.clname, '3C48', band='Q',obsdate='2002-04-20', hosts=hosts
            )
        )
        self.cl.open(self.clname)
        self.assertEqual(self.cl.length(), 385, 'Incorrect number of components')
        ws = self.cl.getkeyword('web_service')
        self.assertEqual(ws['band'], 'Q', 'Incorrect band in web_service metadata')
        self.assertEqual(ws['source'], '3C48', 'Incorrect source in web_service metadata')

       
    def test_direction(self):
        """"Test direction input""""
        hosts = [self.hostname]
        direction = 'J2000 01:37:41.1 33.09.32'
        self.query_server(
            lambda: calmod(
                self.clname, direction=direction, band='Q',obsdate=50000, hosts=hosts
            )
        )
        self.cl.open(self.clname)
        self.assertEqual(self.cl.length(), 385, 'Incorrect number of components')
        ws = self.cl.getkeyword('web_service')
        self.assertEqual(ws['band'], 'Q', 'Incorrect band in web_service metadata')
    """

if __name__ == '__main__':
     unittest.main()
