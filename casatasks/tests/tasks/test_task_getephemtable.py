#########################################################################
# test_task_getephemtable.py
# Copyright (C) 2023
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
# https://casadocs.readthedocs.io/en/stable/api/tt/casatasks.data.getephemtable.html
#
##########################################################################
import contextlib
import csv
import os
import shutil
import tempfile
import unittest
from unittest.mock import Mock
from unittest.mock import patch
from urllib.error import URLError, HTTPError
import certifi
import pytest
import numpy as np

from casatestutils import testhelper as th

from casatasks import getephemtable 
from casatasks.private import jplhorizons_query
from casatools import ctsys, table, measures

_tb = table()
_me = measures()

datapath = ctsys.resolve('/unittest/getephemtable/')

class getephemtable_test(unittest.TestCase):
    def setUp(self):
        self.hostname = 'https://ssd.jpl.nasa.gov/api/horizons.api'
        self.outfile = 'testephem.tab'
        self.caltimerange = '2023/09/01/20:00~2023/09/04/20:00'
        self.jdtimerange = 'JD2460189.33333~2460189.88542'
        self.mjdtimerange = 'MJD60188.83333~60189.38542'

    def tearDown(self):
        if os.path.exists(self.outfile):
            shutil.rmtree(self.outfile) 

    def isDatabaseURLreachable():
        from urllib.request import urlopen
        from urllib.error import URLError
        import certifi
        import ssl
        context = ssl.create_default_context(cafile=certifi.where())
        try: 
            urlaccess = urlopen(self.hostname,context=context, timeout=60.0)
            if urlacess.getcode() == 200:
                return True
            else:
                return False
        except:
            return False

    def test_invalid_inputs(self):
        """Test task inputs"""
        with self.assertRaisesRegex(ValueError, r'objectname must be specified'):
            getephemtable()

        with self.assertRaisesRegex(ValueError, r'timerange must be specified'):
            getephemtable(objectname='Titan')

        with self.assertRaisesRegex(ValueError, r'outfile must be specified'):
            getephemtable(objectname='Titan', timerange=self.caltimerange)

        with self.assertRaisesRegex(RuntimeError, r'objectname is given as an ID'):
            getephemtable(objectname='606',timerange=self.caltimerange, outfile=self.outfile) 

        with self.assertRaisesRegex(ValueError, r'timerange needs to be specified with'):
            getephemtable(objectname='Titan', timerange='2023/09/01/20:00  2023/09/04/20:00', outfile=self.outfile)

        with self.assertRaisesRegex(ValueError, r'Error translating stop time of timerange specified in JD'):
            getephemtable(objectname='Titan', timerange='JD 2460189.8~ 2023/09/04/00', outfile=self.outfile)

        with self.assertRaisesRegex(ValueError, r'Error translating stop time of timerange specified in MJD'):
            getephemtable(objectname='Titan', timerange='MJD 60189.3~ 2023/09/04/00', outfile=self.outfile)

        with self.assertRaisesRegex(ValueError, r'Error in timerange format'):
            getephemtable(objectname='Titan', timerange='09-01-2023 20:00~ 09-02-2023 09:20', outfile=self.outfile)

        with self.assertRaisesRegex(ValueError, r'must contains integer value and unit'):
            getephemtable(objectname='Titan', timerange=self.mjdtimerange, interval='15', outfile=self.outfile)

    @patch('casatasks.private.jplhorizons_query.queryhorizons')
    def test_webservice_errors(self, mock_query):
        with self.assertRaisesRegex(HTTPError, r'Not Found'):
            mock_query.side_effect = HTTPError(url='127.0.0.1', code=404, hdrs={}, fp=None, msg='Not Found')
            getephemtable(objectname='Titan', timerange=self.caltimerange, outfile=self.outfile)

        with self.assertRaisesRegex(URLError, r'Unknown host'):
            mock_query.side_effect = URLError('Unknown host')
            getephemtable(objectname='Titan', timerange=self.caltimerange, outfile=self.outfile)



    @unittest.skipIf(isDatabaseURLreachable(), "JPL-Horizons data server is not reachable")
    def test_tablegeneration(self):
        """Test ephem table generation"""
        getephemtable(objectname='Titan', timerange=self.caltimerange, outfile=self.outfile)

        self.assertTrue(os.path.exists(self.outfile))
        # make sure the table is readable by measures
        self.assertTrue(_me.framecomet(self.outfile))

    @unittest.skipIf(isDatabaseURLreachable(), "JPL-Horizons data server is not reachable")
    def test_saverawdata(self):
        """Test raw query result saving"""
        getephemtable(objectname='Titan', timerange=self.caltimerange, outfile=self.outfile, rawdatafile='saved_rawqueryresult.txt')
        self.assertTrue(os.path.exists(self.outfile))
        # make sure the table is readable by measures
        self.assertTrue(_me.framecomet(self.outfile))



if __name__ == '__main__':
    unittest.main()
