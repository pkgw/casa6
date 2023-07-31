##########################################################################
# test_tool_image_fitsheader.py
#
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
# Based on the requirements listed in casadocs found here:
# https://casadocs.readthedocs.io/en/latest/api/tt/casatools.image.html
#
# Methods tested in this script: fitsheader
#
##########################################################################

import numpy as np
import unittest
from casatools import image as imtool
from numpy import array, array_equal, allclose, ndarray


def isnumarray( ary ):
    ### check if array is numeric
    return np.issubdtype(ary.dtype, np.number)

class ImageBase(unittest.TestCase):
    def setUp(self):
        self._myim = imtool( )
    def compare_dicts( self, dict1, dict2, exclude_keys=[ ] ):
        diff = { }

        for k, v in dict1.items( ):
            ###
            ###  more conditions would be needed to compare nested dictionaries
            ###
            if k in exclude_keys:
                ### skip excluded keys
                continue
            if isinstance(dict1[k], ndarray) and isinstance(dict2[k], ndarray):
                ### compare arrays
                if isnumarray(dict1[k]) and isnumarray(dict2[k]):
                    # numeric array compare
                    if not allclose(dict1[k],dict2[k]):
                        diff[k] = ( dict1[k], dict2[k] )
                elif not isnumarray(dict1[k]) and not isnumarray(dict2[k]):
                    # string array compare
                    if not array_equal(dict1[k],dict2[k]):
                        diff[k] = ( dict1[k], dict2[k] )
                else:
                    # different sorts of arrays
                    diff[k] = ( dict1[k], dict2[k] )
                continue
            if isinstance(dict1[k], ndarray) or isinstance(dict2[k], ndarray):
                ### one array and one non-array is a difference
                diff[k] = ( dict1[k], dict2[k] )
                continue
            if dict1[k] != dict2[k]:
                ### compare singleton values
                diff[k] = ( dict1[k], dict2[k] )

        return diff

class ia_findsources_test(ImageBase):

    def tearDown(self):
        self._myim.done()

    def test_fitsheader(self):
        """test various units are allowed"""
        self._myim.maketestimage()
        expected_header = { 'BITPIX': -32, 'BMAJ': 0.014861112, 'BMIN': 0.0094999997, 'BPA': 6.0, 'BSCALE': 1.0,
                            'BTYPE': 'Intensity', 'BUNIT': 'Jy/beam ', 'BZERO': 0.0,
                            'CDELT': array([-0.00222222,  0.00333333]), 'CRPIX': array([56., 38.]),
                            'CRVAL': array([0., 0.]), 'CTYPE': array(['RA---SIN', 'DEC--SIN'], dtype='<U8'),
                            'CUNIT1': 'deg     ', 'CUNIT2': 'deg     ', 'DATE': '2023-07-31T17:55:28.349375',
                            'END': '', 'EQUINOX': 2000.0, 'EXTEND': True,
                            'HISTORY': array([ "  File modified by user 'dbarnes' with fv  on 1999-07-28T14:19:41",
                                               "  File modified by user 'dbarnes' with fv  on 1999-07-28T14:21:43",
                                               'CASA START LOGTABLE',
                                               "2023-07-31T17:55:28 INFO SRCCODE='::image::maketestimage'",
                                               'Ran ia.maketestimage',
                                               "2023-07-31T17:55:28 INFO SRCCODE='::image::maketestimage'",
                                               'ia.maketestimage(outfile="", overwrite=false)',
                                               'CASA END LOGTABLE'], dtype='<U65'),
                            'IMAGENME': 'Temporary_Image', 'LATPOLE': 0.0, 'LONPOLE': 180.0,
                            'NAXIS': array([  2, 113,  76]), 'OBJECT': '        ',
                            'ORIGIN': 'casacore-@PROJECT_VERSION@', 'PC1_1': 1.0, 'PC1_2': 0.0, 'PC2_1': -0.0,
                            'PC2_2': 1.0, 'PV2_1': 0.0, 'PV2_2': 0.0, 'RADESYS': 'FK5     ',
                            'SIMPLE': True, 'TIMESYS': 'UTC     ' }

        difference = self.compare_dicts(expected_header, self._myim.fitsheader( ),['HISTORY','DATE'])
        self.assertEqual( len(difference), 0, f'''expected no differences between reference and created fits headers, but found: {repr(difference)}''' )

if __name__ == '__main__':
    unittest.main()

