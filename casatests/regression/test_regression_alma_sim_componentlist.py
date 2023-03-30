########################################################################
# test_regression_alma_sim_componentlist.py
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
# Based on CASAguides
# https://casaguides.nrao.edu/index.php/Simulation_Guide_Component_Lists_(CASA_5.4)
#
##########################################################################

import os
import shutil
import unittest

# On purpose this script does not import any casatool or casatask
# It is meant to run with casa -c testscript using a CASA monolithic tarball,
# which will verify if tasks are automatically imported in the tarball

class ALMASimComponentListTest(unittest.TestCase):

    def setUp(self) -> None:
        self.imname = 'Gaussian.im'
        self.imfits = 'Gaussian.fits'
        self.point = 'point.cl'
        self.imagelist = 'Image_list'

    def tearDown(self) -> None:
        # the .png files created by simobserve and simanalyse will be kept
        # inside a directory called png_list. These png images will be available
        # at the Artifacts tab of the Bamboo plan that runs these tests
        shutil.rmtree(self.imname)
        os.remove(self.imfits)
        shutil.rmtree(self.point)
        if os.path.exists('png_list'):
            shutil.rmtree('png_list')
        os.system('mkdir png_list')
        os.system('cp '+self.imagelist+'/*.png png_list')
        shutil.rmtree(self.imagelist)
        cl.done()
        ia.done()

    def test_simulation_with_componentlist(self):
        """Make FITS image and simulate observations with component list"""
        direction = "J2000 10h00m00.0s -30d00m00.0s"
        cl.done()
        cl.addcomponent(dir=direction, flux=1.0, fluxunit='Jy', freq='230.0GHz', shape="Gaussian",
                        majoraxis="0.1arcmin", minoraxis='0.05arcmin', positionangle='45.0deg')

        ia.fromshape(self.imname,[256,256,1,1],overwrite=True)
        self.assertTrue(os.path.exists(self.imname),"Cannot find %s"%self.imname)

        cs = ia.coordsys()
        cs.setunits(['rad','rad','','Hz'])
        cell_rad = qa.convert(qa.quantity("0.1arcsec"),"rad")['value']
        cs.setincrement([-cell_rad,cell_rad],'direction')
        cs.setreferencevalue([qa.convert("10h",'rad')['value'],qa.convert("-30deg",'rad')['value']],type="direction")
        cs.setreferencevalue("230GHz",'spectral')
        cs.setincrement('1GHz','spectral')
        ia.setcoordsys(cs.torecord())
        ia.setbrightnessunit("Jy/pixel")
        ia.modify(cl.torecord(),subtract=False)
        ia.done()

        # Export the image to a FITS file
        exportfits(imagename=self.imname,fitsimage=self.imfits,overwrite=True)
        self.assertTrue(os.path.exists(self.imname),"Cannot find %s"%self.imname)

        cl.done()
        cl.addcomponent(dir="J2000 10h00m00.08s -30d00m02.0s", flux=0.1, fluxunit='Jy', freq='230.0GHz', shape="point")
        cl.addcomponent(dir="J2000 09h59m59.92s -29d59m58.0s", flux=0.1, fluxunit='Jy', freq='230.0GHz', shape="point")
        cl.addcomponent(dir="J2000 10h00m00.40s -29d59m55.0s", flux=0.1, fluxunit='Jy', freq='230.0GHz', shape="point")
        cl.addcomponent(dir="J2000 09h59m59.60s -30d00m05.0s", flux=0.1, fluxunit='Jy', freq='230.0GHz', shape="point")
        cl.rename(self.point)
        cl.done()
        self.assertTrue(os.path.exists(self.point),"Cannot find %s"%self.point)

        # Simulate observations of the point sources the Gaussian flux distribution in the FITS file
        simobserve(project = self.imagelist, skymodel = self.imfits, inwidth = "1GHz", complist = self.point,
                   compwidth = '1GHz', direction = "J2000 10h00m00.0s -30d00m00.0s", obsmode = "int",
                   antennalist = 'alma.cycle6.1.cfg', totaltime = "28800s", thermalnoise = '')

        self.assertTrue(os.path.exists(self.imagelist),"Cannot find %s"%self.imagelist)
        self.assertTrue(os.path.exists(os.path.join(self.imagelist,"Image_list.alma.cycle6.1.ms")))

        # Make a map of the emision using the project created in the previous simobser step
        simanalyze(project = self.imagelist, imsize = [256,256], imdirection = "J2000 10h00m00.0s -30d00m00.0s",
                   cell = '0.1arcsec', niter = 5000, threshold = '10.0mJy/beam', analyze = True)

        self.assertTrue(os.path.exists(os.path.join(self.imagelist,"Image_list.alma.cycle6.1.image")))

if __name__ == "__main__":
    unittest.main()

