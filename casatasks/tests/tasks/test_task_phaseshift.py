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

from casatools import (
        componentlist, ctsys, image, measures, ms, msmetadata, quanta,
        regionmanager, simulator, table
    )
from casatasks import flagdata, imstat, phaseshift, tclean
from casatasks.private import simutil
cl = componentlist()
ia = image()
md = msmetadata()
me = measures()
ms = ms()
qa = quanta()
rg = regionmanager()
sm = simulator()
tb = table()
datadir = os.path.join('unittest', 'phaseshift')
ctsys_resolve = ctsys.resolve

datadir = os.path.join('unittest', 'phaseshift')
datapath = ctsys_resolve(os.path.join(datadir, 'refim_twopoints_twochan.ms'))
datapath_Itziar = ctsys_resolve(os.path.join(datadir, 'Itziar.ms'))
datapath_ngc = ctsys_resolve(os.path.join(datadir, 'ngc7538_ut.ms'))
datapath_nep = ctsys_resolve(os.path.join(datadir, 'nep2-shrunk.ms'))
datapath_mms = ctsys_resolve(
    os.path.join(datadir, 'uid___X02_X3d737_X1_01_small.mms')
)


def change_perms(path):
    os.chmod(path, 0o777)
    for root, dirs, files in os.walk(path):
        for d in dirs:
            os.chmod(os.path.join(root, d), 0o777)
        for f in files:
            os.chmod(os.path.join(root, f), 0o777)


datacopy = 'datacopy.ms'
datacopy_Itziar = 'Itziar_copy.ms'
datacopy_ngc = 'ngc_copy.ms'
datacopy_nep = 'nep_copy.ms'
datacopy_mms = 'mms_copy.mms'
output = 'phaseshiftout.ms'


class phaseshift_base_checks(unittest.TestCase):

    def check_nrows(self, vis, expected_nrows):
        """ Simple check to ensure the correct number of rows has been produced in the main
        table of an (output) MS. Checks the vis exists on disk and looks at the number of
        rows of the TIME column. """
        tbt = table()
        try:
            tbt.open(vis)
            nrows = len(tbt.getcol('TIME'))
        finally:
            tbt.close()

        self.assertTrue(os.path.exists(vis), f"MeasurementSet not found: {vis}")
        self.assertEqual(nrows, expected_nrows,
                         f"Incorrect number of rows found in MS: {vis}")

    def check_field_subtable(self, outputvis, inputvis, new_centers,
                             field_selection=None):
        """
        Ensures that the field subtable of outputvis has the expected phase centers
        (from 'new_centers' passed to the phasecenter parameter of phaseshift),
        comparing also unmodified ('passthrough') fields with the 'inputvis'.

        :param outputvis: output, phase-shifted MS to check
        :param inputvis: input MS, to compare against when some fields unchanged
        :param new_centers: string or dict with new phase center(s)
        :param field_selection: field selected as integer if used (only single value supported)

        """
        from casatasks.private.task_phaseshift import _convert_to_ra_dec_j2000

        def get_expected_output_ra_dec(row, new_centers, input_phase_col,
                                       field_selection):
            field_id = str(row)
            if isinstance(new_centers, dict):
                if field_id in new_centers:
                    ra_rad, dec_rad = _convert_to_ra_dec_j2000(
                        new_centers[field_id])
                else:
                    ra_rad = input_phase_col[0, 0, row]
                    dec_rad = input_phase_col[1, 0, row]
            else:
                if field_selection and row != field_selection:
                    ra_rad = input_phase_col[0, 0, row]
                    dec_rad = input_phase_col[1, 0, row]
                else:
                    ra_rad, dec_rad = _convert_to_ra_dec_j2000(new_centers)

            return ra_rad, dec_rad

        def get_field_subt_phasedir_col(vis_path):
            try:
                tblocal = table()
                tblocal.open(vis_path + '/FIELD', nomodify=True)
                phase_col = tblocal.getcol('PHASE_DIR')
            finally:
                tblocal.done()

            return phase_col


        phase_col = get_field_subt_phasedir_col(outputvis)
        input_phase_col = get_field_subt_phasedir_col(inputvis)
        for row in range(0, phase_col.shape[-1]):
            ra_rad, dec_rad = get_expected_output_ra_dec(row, new_centers,
                                                         input_phase_col,
                                                         field_selection)

            # The 0 in the middle is the 'NUM_POLY' axis
            self.assertEqual(phase_col[0, 0, row], ra_rad,
                             f"unexpected PHASE_DIR ra value in row {row} "
                             f"(with {new_centers=})")
            self.assertEqual(phase_col[1, 0, row], dec_rad,
                             f"unexpected PHASE_DIR dec value in row {row} "
                             f"(with {new_centers=})")


class phaseshift_test(phaseshift_base_checks):

    def setUp(self):
        shutil.copytree(datapath, datacopy)
        shutil.copytree(datapath_Itziar, datacopy_Itziar)
        shutil.copytree(datapath_ngc, datacopy_ngc)
        shutil.copytree(datapath_nep, datacopy_nep)
        shutil.copytree(datapath_mms, datacopy_mms)

        change_perms(datacopy)
        change_perms(datacopy_Itziar)
        change_perms(datacopy_ngc)
        change_perms(datacopy_nep)
        change_perms(datacopy_mms)

    def tearDown(self):
        shutil.rmtree(datacopy)
        shutil.rmtree(datacopy_Itziar)
        shutil.rmtree(datacopy_ngc)
        shutil.rmtree(datacopy_nep)
        shutil.rmtree(datacopy_mms)

        if os.path.exists('post_phaseshift.ms'):
            shutil.rmtree('post_phaseshift.ms')

        if os.path.exists(output):
            shutil.rmtree(output)

    def test_takesVis(self):
        ''' Check that the task requires a valid input MS '''
        phasecenter = 'J2000 19h53m50 40d06m00'
        result = phaseshift(
            datacopy, outputvis=output,
            phasecenter=phasecenter
        )
        self.check_nrows(output, 63180)
        self.check_field_subtable(output, datacopy, phasecenter)

    def test_outvis(self):
        '''
        Check that the outvis parameter specifies the name of the output
        '''
        phasecenter = 'J2000 19h53m50 40d06m00'
        phaseshift(
            datacopy, outputvis=output,
            phasecenter=phasecenter
        )

        self.check_nrows(output, 63180)
        self.check_field_subtable(output, datacopy, phasecenter)

    def test_fieldSelect(self):
        ''' Check the field selection parameter '''
        phasecenter = 'J2000 00h00m01 -29d55m40'
        phaseshift(
            datacopy_Itziar, outputvis=output,
            phasecenter=phasecenter, field='2'
        )

        self.check_nrows(output, 6125)
        self.check_field_subtable(output, datacopy_Itziar, phasecenter,
                                  field_selection=2)

    def test_spwSelect(self):
        ''' Check the spw selection parameter '''
        phasecenter = 'B1950_VLA 23h11m54 61d10m54'
        phaseshift(
            datacopy_ngc, outputvis=output,
            phasecenter=phasecenter, spw='1'
        )
        tb.open(output)
        data_selected = len(tb.getcol('TIME'))
        tb.close()

        self.check_nrows(output, 13338)
        self.check_field_subtable(output, datacopy_ngc, phasecenter)

    def test_intentSelect(self):
        ''' Check the intent selection parameter '''
        phasecenter = 'ICRS 00h06m14 -06d23m35'
        phaseshift(
            datacopy_nep, outputvis=output,
            phasecenter=phasecenter, intent='*FLUX*'
        )

        self.check_nrows(output, 570)
        self.check_field_subtable(output, datacopy_nep, phasecenter)

    def test_arraySelect(self):
        ''' Check the array selection parameter '''
        msg = "specified array incorrectly found"
        with self.assertRaises(RuntimeError, msg=msg):
            phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter='ICRS 00h06m14 -06d23m35',
                    array='1'
            )

        phasecenter = 'ICRS 00h06m14 -06d23m35'
        phaseshift(
            datacopy_nep, outputvis=output,
            phasecenter=phasecenter, array='0'
        )

        self.check_nrows(output, 6270)
        self.check_field_subtable(output, datacopy_nep, phasecenter)

    def test_observationSelect(self):
        ''' Check the observation selection parameter '''
        msg = "Observation not out of range"
        with self.assertRaises(RuntimeError, msg=msg):
            phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter='ICRS 00h06m14 -06d23m35', observation='1'
            )

        phasecenter = 'ICRS 00h06m14 -06d23m35'
        phaseshift(
            datacopy_nep, outputvis=output,
            phasecenter=phasecenter, observation='0'
        )

        self.check_nrows(output, 6270)
        self.check_field_subtable(output, datacopy_nep, phasecenter)

    def test_keepsMMS(self):
        '''
        Test the keepmms paramter creates the output as an MMS
        if the input is one as well
        '''
        phasecenter = 'J2000 05h30m48 13d31m48'
        phaseshift(
            datacopy_mms, outputvis=output,
            phasecenter=phasecenter, keepmms=False
        )
        ms.open(output)
        is_mms = ms.ismultims()
        ms.close()

        self.assertFalse(is_mms)
        self.check_nrows(output, 1080)
        self.check_field_subtable(output, datacopy_mms, phasecenter)

    def test_datacolumn(self):
        '''
        Check that this parameter selects which datacolumns to write
        to the output MS
        '''
        msg = "Data column incorrectly present"
        with self.assertRaises(RuntimeError, msg=msg):
            phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter='ICRS 00h06m14 -06d23m35',
                    datacolumn='MODEL'
            )

        # running to completion indicates success in CASA 6
        phasecenter = 'ICRS 00h06m14 -06d23m35'
        phaseshift(
            datacopy_nep, outputvis=output,
            phasecenter=phasecenter, datacolumn='DATA'
        )

        self.check_nrows(output, 6270)
        self.check_field_subtable(output, datacopy_nep, phasecenter)

    def test_phasecenter(self):
        '''
        Check that this parameter sets the sky coordinates of the new
        phasecenter
        '''
        msg = 'Empty phasecenter param incorrectly runs'
        with self.assertRaises(ValueError, msg=msg):
            phaseshift(
                    datacopy_nep, outputvis=output,
                    phasecenter=''
            )

        phasecenter = 'ICRS 00h06m14 -08d23m35'
        phaseshift(datacopy_nep, outputvis=output, phasecenter=phasecenter)

        self.check_nrows(output, 6270)
        tb.open(output)
        data_mean = np.mean(tb.getcol('DATA'))
        tb.close()

        self.assertTrue(np.isclose(
            data_mean, -0.00968202886279957-0.004072808512879953j)
        )
        self.check_field_subtable(output, datacopy_nep, phasecenter)

    def test_shiftAndCompare(self):
        '''
        Check that changing the phasecenter with phaseshift and
        reverting with tclean results in the correct flux values at
        selected pixel locations
        '''
        # Run phaseshift to shift the MS phasecenter to a new location.
        post_vis = 'post_phaseshift.ms'
        os.system('rm -rf ' + post_vis)
        phasecenter = 'J2000 19h53m50 40d06m00'
        phaseshift(
            vis=datacopy, outputvis=post_vis,
            phasecenter=phasecenter
        )
        self.check_nrows(post_vis, 63180)
        self.check_field_subtable(post_vis, datacopy, phasecenter)

        # (1) Imaging on the original dataset
        os.system('rm -rf im2_pre*')
        tclean(
            vis=datacopy, imagename='im2_pre',
            imsize=2048, cell='5arcsec', niter=0,
            gridder='wproject', wprojplanes=128, pblimit=-0.1,
            phasecenter='J2000 19h59m28.449 40d44m01.199'
        )

        # (2) Image the phaseshifted dataset at it's new phasecenter as
        # image center.
        post_image = ''
        post_image = 'im2_post_phaseshift'
        os.system('rm -rf im2_post_phaseshift*')

        tclean(
            vis=post_vis, imagename='im2_post_phaseshift',
            imsize=2048, cell='5arcsec', niter=0,
            gridder='wproject', wprojplanes=128, pblimit=-0.1
        )

        # (3) Imaging on phaseshifted dataset, but with the imaging phasecenter
        # set. If this is working correctly, it should shift back to the same
        # source positions as the previous tclean result.
        post_image = 'im2_post_phaseshift_tclean_phasecenter'
        os.system('rm -rf im2_post_phaseshift_tclean_phasecenter*')

        tclean(
            vis=post_vis, imagename=post_image,
            imsize=2048, cell='5arcsec', niter=0, gridder='wproject',
            wprojplanes=128, pblimit=-0.1,
            phasecenter='J2000 19h59m28.449 40d44m01.199'
        )

        # In the above 3 images, (1) has the correct locations.
        # Both (2) and (3) show the offset error when viewed in world
        # coordinates. Open in the viewer as an image stack, and step
        # through.
        # For comparisons with (1), pick the result from (3) because when
        # this works correctly the sources should appear at the same
        # pixel location as in (1).  This test is encoded below.

        ia.open('im2_pre.image')
        src1_pre = ia.pixelvalue([1024, 1024, 0, 0])['value']
        src2_pre = ia.pixelvalue([1132, 1168, 0, 0])['value']
        ia.close()
        ia.open(post_image+'.image')
        src1_post = ia.pixelvalue([1024, 1024, 0, 0])['value']
        src2_post = ia.pixelvalue([1132, 1168, 0, 0])['value']
        ia.close()

        print("Image value at source locations")
        print("Original MS : "+str(src1_pre) + " and " + str(src2_pre))
        print("Phase shifted MS : "+str(src1_post) + " and " + str(src2_post))

        os.system('rm -rf im2_pre*')
        os.system('rm -rf im2_post_phaseshift*')
        os.system('rm -rf im2_post_phaseshift_tclean_phasecenter*')

        self.assertTrue(
            np.isclose(src1_pre['value'], src1_post['value'], rtol=0.01)
        )
        self.assertTrue(
            np.isclose(src2_pre['value'], src2_post['value'], rtol=0.01)
        )


class reference_frame_tests(unittest.TestCase):
    # much of this code is adapted from that of rurvashi
    # https://gitlab.nrao.edu/rurvashi/simulation-in-casa-6/-/blob/master/Simulation_Script_Demo.ipynb

    comp_list = 'sim_onepoint.cl'
    orig_ms = 'sim_data.ms'
    orig_im = 'im_orig'
    pshift_ms = 'sim_data_pshift.ms'
    pshift_im = 'im_post_phaseshift'
    pshift_shiftback_im = 'im_post_phaseshift_tclean_shiftback'
    exp_flux = 5

    @classmethod
    def __delete_intermediate_products(cls):
        if os.path.exists(cls.pshift_ms):
            shutil.rmtree(cls.pshift_ms)
        for im in (cls.pshift_im, cls.pshift_shiftback_im):
            for path in glob.glob(im + '*'):
                shutil.rmtree(path)

    @classmethod
    def __delete_data(cls):
        cls.__delete_intermediate_products()
        for x in (cls.orig_ms, cls.comp_list):
            if os.path.exists(x):
                shutil.rmtree(x)
        for path in glob.glob(cls.orig_im + '*'):
            shutil.rmtree(path)

    def setUp(self):
        self.__delete_data()

    def tearDown(self):
        self.__delete_data()

    @classmethod
    def tearDownClass(cls):
        cls.__delete_data()

    @classmethod
    def __phase_center_string(cls, ra, dec, frame):
        return ' '.join([frame, ra, dec])

    @classmethod
    def __makeMSFrame(cls, radir, decdir, dirframe):
        """
        Construct an empty Measurement Set that has the desired
        observation setup.
        """
        # Open the simulator
        sm.open(ms=cls.orig_ms)

        # Read/create an antenna configuration.
        # Canned antenna config text files are located at
        # /home/casa/data/trunk/alma/simmos/*cfg
        antennalist = os.path.join(
            ctsys.resolve("alma/simmos"), "vla.d.cfg"
        )


        # Fictitious telescopes can be simulated by specifying x, y, z, d,
        # an, telname, antpos.
        # x,y,z are locations in meters in ITRF (Earth centered)
        # coordinates.
        # d, an are lists of antenna diameter and name.
        # telname and obspos are the name and coordinates of the
        # observatory.
        (x, y, z, d, an, an2, telname, obspos) = (
            simutil.simutil().readantenna(antennalist)
        )

        # Set the antenna configuration
        sm.setconfig(
            telescopename=telname, x=x, y=y, z=z, dishdiameter=d,
            mount=['alt-az'], antname=an, coordsystem='global',
            referencelocation=me.observatory(telname)
        )

        # Set the polarization mode (this goes to the FEED subtable)
        sm.setfeed(mode='perfect R L', pol=[''])

        # Set the spectral window and polarization (one
        # data-description-id).
        # Call multiple times with different names for multiple SPWs or
        # pol setups.
        sm.setspwindow(
            spwname="LBand", freq='1.0GHz', deltafreq='0.1GHz',
            freqresolution='0.2GHz', nchannels=1, stokes='RR LL'
        )

        # Setup source/field information (i.e. where the observation phase
        # center is) Call multiple times for different pointings or source
        # locations.
        sm.setfield(
            sourcename="fake", sourcedirection=me.direction(
                rf=dirframe, v0=radir, v1=decdir
            )
        )

        # Set shadow/elevation limits (if you care). These set flags.
        sm.setlimits(shadowlimit=0.01, elevationlimit='1deg')

        # Leave autocorrelations out of the MS.
        sm.setauto(autocorrwt=0.0)

        # Set the integration time, and the convention to use for timerange
        # specification
        # Note : It is convenient to pick the hourangle mode as all times
        #   specified in sm.observe() will be relative to when the source
        #   transits.
        sm.settimes(
            integrationtime='2000s', usehourangle=True,
            referencetime=me.epoch('UTC', '2019/10/4/00:00:00')
        )

        # Construct MS metadata and UVW values for one scan and ddid
        # Call multiple times for multiple scans.
        # Call this with different sourcenames (fields) and spw/pol
        # settings as defined above.
        # Timesteps will be defined in intervals of 'integrationtime',
        # between starttime and stoptime.
        sm.observe(
            sourcename="fake", spwname='LBand', starttime='-5.0h',
            stoptime='+5.0h'
        )
        # Close the simulator
        sm.close()
        # Unflag everything (unless you care about elevation/shadow flags)
        flagdata(vis=cls.orig_ms, mode='unflag')

    @classmethod
    def __makeCompList(cls, ra, dec, frame):
        # Add sources, one at a time.
        # Call multiple times to add multiple sources.
        # ( Change the 'dir', obviously )
        cl.addcomponent(
            dir=cls.__phase_center_string(ra, dec, frame),
            flux=cls.exp_flux,      # For a gaussian, this is the
                                    # integrated area.
            fluxunit='Jy', freq='1.5GHz', shape='point',
            spectrumtype="constant"
        )
        # Save the file
        cl.rename(filename=cls.comp_list)
        cl.done()

    @classmethod
    def __sim2fields(cls, radir, decdir, dirframe, offset):
        """
        Construct an empty Measurement Set with two fieldsthat has
        the desired observation setup.
        """
        # Open the simulator
        sm.open(ms=cls.orig_ms)

        # Read/create an antenna configuration.
        # Canned antenna config text files are located at
        # /home/casa/data/trunk/alma/simmos/*cfg
        antennalist = os.path.join(
            ctsys.resolve("alma/simmos"), "vla.d.cfg"
        )


        # Fictitious telescopes can be simulated by specifying x, y, z, d,
        # an, telname, antpos.
        # x,y,z are locations in meters in ITRF (Earth centered)
        # coordinates.
        # d, an are lists of antenna diameter and name.
        # telname and obspos are the name and coordinates of the
        # observatory.
        (x, y, z, d, an, an2, telname, obspos) = (
            simutil.simutil().readantenna(antennalist)
        )

        # Set the antenna configuration
        sm.setconfig(
            telescopename=telname, x=x, y=y, z=z, dishdiameter=d,
            mount=['alt-az'], antname=an, coordsystem='global',
            referencelocation=me.observatory(telname)
        )

        # Set the polarization mode (this goes to the FEED subtable)
        sm.setfeed(mode='perfect R L', pol=[''])

        # Set the spectral window and polarization (one
        # data-description-id).
        # Call multiple times with different names for multiple SPWs or
        # pol setups.
        sm.setspwindow(
            spwname="LBand", freq='1.0GHz', deltafreq='0.1GHz',
            freqresolution='0.2GHz', nchannels=1, stokes='RR LL'
        )

        # Setup source/field information (i.e. where the observation phase
        # center is) Call multiple times for different pointings or source
        # locations.
        sm.setfield(
            sourcename="fake", sourcedirection=me.direction(
                rf=dirframe, v0=radir, v1=decdir
            )
        )
        # the second field is 10deg north of the first, so should be
        # emission free
        sm.setfield(
            sourcename="pretend", sourcedirection=me.direction(
                rf=dirframe, v0=radir, v1=qa.tos(
                    qa.add(qa.quantity(decdir), qa.quantity(offset))
                )
            )
        )
        # Set shadow/elevation limits (if you care). These set flags.
        sm.setlimits(shadowlimit=0.01, elevationlimit='1deg')

        # Leave autocorrelations out of the MS.
        sm.setauto(autocorrwt=0.0)

        # Set the integration time, and the convention to use for timerange
        # specification
        # Note : It is convenient to pick the hourangle mode as all times
        #   specified in sm.observe() will be relative to when the source
        #   transits.
        sm.settimes(
            integrationtime='2000s', usehourangle=True,
            referencetime=me.epoch('UTC', '2019/10/4/00:00:00')
        )

        # Construct MS metadata and UVW values for one scan and ddid
        # Call multiple times for multiple scans.
        # Call this with different sourcenames (fields) and spw/pol
        # settings as defined above.
        # Timesteps will be defined in intervals of 'integrationtime',
        # between starttime and stoptime.
        sm.observe(
            sourcename="fake", spwname='LBand', starttime='-5.0h',
            stoptime='-2.5h'
        )
        sm.observe(
            sourcename="pretend", spwname='LBand', starttime='-2.5h',
            stoptime='0h'
        )
        sm.observe(
            sourcename="fake", spwname='LBand', starttime='0h',
            stoptime='2.5h'
        )
        sm.observe(
            sourcename="pretend", spwname='LBand', starttime='2.5h',
            stoptime='5h'
        )

        # Close the simulator
        sm.close()
        # Unflag everything (unless you care about elevation/shadow flags)
        flagdata(vis=cls.orig_ms, mode='unflag')

    @classmethod
    def __predictSimFromComplist(cls):
        sm.openfromms(cls.orig_ms)
        # Predict from a component list
        sm.predict(complist=cls.comp_list, incremental=False)
        sm.close()

    @classmethod
    def __createImage(cls, msname, imagename, phasecenter):
        for path in glob.glob(imagename + '*'):
            shutil.rmtree(path)
        tclean(
            vis=msname, imagename=imagename, datacolumn='data',
            imsize=256, cell='8.0arcsec', gridder='standard',
            niter=20, gain=0.3, pblimit=-0.1,
            phasecenter=phasecenter
        )

    def __compare(self, imagename, radir, decdir, dirframe):
        ia.open('.'.join([imagename, 'image']))
        stats = ia.statistics()
        maxpos = stats['maxpos']
        (xc, yc) = maxpos[0:2]
        blc = [xc-10, yc-10, 0, 0]
        trc = [xc+10, yc+10, 0, 0]
        fit = ia.fitcomponents(region=rg.box(blc=blc, trc=trc))
        ia.done()
        cl.fromrecord(fit['deconvolved'])
        pos = cl.getrefdir(0)
        flux = cl.getfluxvalue(0)
        cl.done()
        expec = me.direction(dirframe, radir, decdir)
        diff = me.separation(pos, expec)
        self.assertTrue(
            qa.lt(diff, qa.quantity('0.15arcsec')),
            'position difference is too large for ' + str(pos)
            + ': ' + qa.tos(qa.convert(diff, 'arcsec'))
        )
        self.assertAlmostEqual(
            flux[0], self.exp_flux,
            msg='flux differs by too much: got: ' + str(flux[0])
            + ' expected: ' + str(self.exp_flux), delta=0.01
        )

    def __compare_ms(self, dirframe):
        tb.open(
            ctsys_resolve(
                os.path.join(
                    datadir, 'phaseshift_test_frames_expected.ms'
                )
            )
        )
        expuvw = tb.getcol('UVW')
        tb.done()
        tb.open(self.pshift_ms)
        gotuvw = tb.getcol('UVW')
        gotdata = tb.getcol('DATA')
        tb.done()
        # absolute tolerance in meters
        self.assertTrue(
            np.allclose(gotuvw, expuvw, atol=3e-6),
            'UVW do not match for ' + dirframe
        )
        # point source at phase center, so data should have
        # const amplitude and zero phase
        self.assertTrue(
            np.allclose(np.abs(gotdata), self.exp_flux, rtol=1e-6),
            'Amplitudes do not match for ' + dirframe
        )
        # phases in radians
        self.assertTrue(
            np.allclose(np.angle(gotdata), 0, atol=1e-4),
            'Phases do not match for ' + dirframe
        )

    def __run_direction_test(self, p, radir, decdir, dirframe):
        pr = [p['lon'], qa.time(p['lon'], 10)[0]]
        pd = [p['lat'], qa.time(p['lat'], 10)[0]]
        for unit in ['deg', 'rad']:
            pr.append(qa.tos(qa.convert(qa.toangle(p['lon']), unit)))
            pd.append(qa.tos(qa.convert(qa.toangle(p['lat']), unit)))
        for lon in pr:
            for lat in pd:
                shifted_pcenter = self.__phase_center_string(
                    lon, lat, p['frame']
                )
                # run phaseshift
                phaseshift(
                    vis=self.orig_ms, outputvis=self.pshift_ms,
                    phasecenter=shifted_pcenter
                )
                self.__compare_ms(p['frame'])
                # create image from phaseshifted MS
                self.__createImage(self.pshift_ms, self.pshift_im, "")
                self.__compare(self.pshift_im, radir, decdir, dirframe)
                x = imstat(self.pshift_im + '.image')
                # source should be at image center after phase shift
                self.assertTrue(
                    ((x['maxpos'] == [128, 128, 0, 0]).all()),
                    msg='Source stats' + str(x) + ' for ' + str(p)
                )
                self.__delete_intermediate_products()

    def test_frames(self):
        # This is the source position
        radir = '19h53m50'
        decdir = '40d06m00'
        dirframe = 'J2000'
        # this is the field center
        fra = '19h59m28.5'
        fdec = '+40.40.01.5'
        fframe = 'J2000'

        def create_input_ms():
            # this is how self.orig_ms was created. Do not delete
            # this or related code even if the data set is now
            # stored in the data repos. The code is useful to have
            # as a record.
            # make the MS
            self.__makeMSFrame(fra, fdec, fframe)
            # Make the component list
            self.__makeCompList(radir, decdir, dirframe)
            # Predict Visibilities
            self.__predictSimFromComplist()

        # create_input_ms()
        shutil.copytree(
            ctsys_resolve(
                os.path.join(datadir, 'phaseshift_test_frames_input.ms')
            ), self.orig_ms
        )
        # image simulated MS, the source is offset from the phase center
        # of the image
        tclean(
            vis=self.orig_ms, imagename=self.orig_im, datacolumn='data',
            imsize=2048, cell='8.0arcsec', gridder='wproject',
            niter=20, gain=0.3, wprojplanes=128, pblimit=-0.1
        )
        x = imstat(self.orig_im + '.image')
        self.assertTrue((x['maxpos'] == [1509, 773, 0, 0]).all())
        # self.__createImage(self.orig_ms, self.orig_im, orig_pcenter)
        self.__compare(self.orig_im, radir, decdir, dirframe)
        # J2000
        j2000 = {'lon': '19h53m50', 'lat': '40d06m00', 'frame': 'J2000'}
        # ICRS coordinates of the above
        icrs = {'lon': '19h53m49.9980', 'lat': '40d06m0.0019', 'frame': 'ICRS'}
        # GALACTIC coordinates of the above
        galactic = {
            'lon': '05h00m21.5326', 'lat': '+006d21m09.7433',
            'frame': 'GALACTIC'
        }
        # B1950_VLA coordinates of the above
        b1950_vla = {
            'lon': '19h52m05.65239', 'lat': '39d58m05.8512',
            'frame': 'B1950_VLA'
        }
        for p in (j2000, icrs, galactic, b1950_vla):
            self.__run_direction_test(p, radir, decdir, dirframe)
        self.__delete_data()

    def test_field(self):
        """Test that a field is correctly chosen in a multi-field MS"""
        def shift_and_clean(myfield, expdir):
            phaseshift(
                vis=self.orig_ms, outputvis=self.pshift_ms,
                phasecenter=pcenter, field=myfield
            )

            field_id = ''
            if len(myfield) > 0:
                try:
                    field_id = int(myfield)
                except ValueError:
                    md.open(self.orig_ms)
                    field_id = md.fieldsforname(myfield)[0]
                    md.done()
            if field_id:
                separation_field = field_id
            else:
                separation_field = 0

            md.open(self.pshift_ms)
            # re-indexing disabled. Output FIELD subtable has all the original fields
            exp_nfields = 2
            self.assertEqual(
                md.nfields(), exp_nfields,
                msg='Wrong number of fields in FIELD subtable for field ' + myfield
            )
            sep = me.separation(md.refdir(field=separation_field), expdir)
            md.done()
            self.assertEqual(
                qa.getvalue(sep), 0,
                msg='Ref direction is wrong for field ' + myfield
                + ' separation is ' + qa.tos(qa.convert(sep, 'arcsec'))
            )

            # check times and baselines
            if field_id == 0:
                exp_ms = ctsys_resolve(
                    os.path.join(
                        datadir, 'phaseshift_test_field_0_expected.ms'
                    )
                )
            elif field_id == 1:
                exp_ms = ctsys_resolve(
                    os.path.join(
                        datadir, 'phaseshift_test_field_1_expected.ms'
                    )
                )
            elif len(field_id) == 0:
                exp_ms = ctsys_resolve(
                    os.path.join(
                        datadir, 'phaseshift_test_field_0_1_expected.ms'
                    )
                )
            tb.open(exp_ms)
            # myfilter = '' if len(myfield) == 0 else 'FIELD_ID='
            # + str(field_id)
            # x = tb.query(myfilter, columns='TIME, ANTENNA1, ANTENNA2')
            exptime = tb.getcol('TIME')
            expant1 = tb.getcol('ANTENNA1')
            expant2 = tb.getcol('ANTENNA2')
            expuvw = tb.getcol('UVW')
            expdata = tb.getcol('DATA')
            # x.done()
            tb.done()
            tb.open(self.pshift_ms)
            gottime = tb.getcol('TIME')
            gotant1 = tb.getcol('ANTENNA1')
            gotant2 = tb.getcol('ANTENNA2')
            gotuvw = tb.getcol('UVW')
            gotdata = tb.getcol('DATA')
            tb.done()
            self.assertTrue(
                (gottime == exptime).all(),
                msg='Failed TIME column test for "' + myfield + '"'
            )
            self.assertTrue(
                (gotant1 == expant1).all(),
                msg='Failed ANTENNA1 column test for "' + myfield + '"'
            )
            self.assertTrue(
                (gotant2 == expant2).all(),
                msg='Failed ANTENNA2 column test for "' + myfield + '"'
            )
            self.assertTrue(
                np.allclose(gotuvw, expuvw),
                msg='Failed UVW column test for "' + myfield + '"'
            )
            self.assertTrue(
                (gotdata == expdata).all(),
                msg='Failed DATA column test for "' + myfield + '"'
            )
            tclean(
                vis=self.pshift_ms, imagename=self.pshift_im,
                datacolumn='data', imsize=256, cell='8.0arcsec',
                gridder='standard', niter=20, gain=0.3, pblimit=-0.1
            )

        # This is the source position
        radir = '19h53m50'
        decdir = '40d06m00'
        dirframe = 'J2000'
        offset = '10deg'

        def create_ms():
            # do not delete this code, even if the MS is now in the data
            # repos; the code is useful to have as a record and a guide.
            # make the MS
            self.__sim2fields(radir, decdir, dirframe, offset)
            # Make the component list
            self.__makeCompList(radir, decdir, dirframe)
            # Predict Visibilities
            self.__predictSimFromComplist()

        # create_ms()
        shutil.copytree(
            ctsys_resolve(
                os.path.join(datadir, 'phaseshift_test_field_input.ms')
            ), self.orig_ms
        )
        # shift first field by 4 pixels north of source
        pcenter = self.__phase_center_string(
            radir,
            qa.tos(qa.add(qa.quantity(decdir), qa.quantity("32arcsec"))),
            dirframe
        )
        expdir = me.direction(dirframe, radir, decdir)
        # test both incarnations of the first field as well as both
        # fields together ('')
        for myfield in ('0', 'fake', ''):
            shift_and_clean(myfield, expdir)
            x = imstat(self.pshift_im + '.image')
            self.assertTrue(
                (x['maxpos'] == [128, 124, 0, 0]).all(),
                msg='maxpos is incorrect'
            )
            self.assertTrue(
                np.isclose(x['max'][0], self.exp_flux, 1e-6),
                msg='max is incorrect, expected ' + str(self.exp_flux)
                + ' got ' + str(x['max'][0])
            )
            self.__delete_intermediate_products()
        # 4 pixel shift of second field, which contains no signal, but
        # just sidelobes of source 10 degrees away
        decref = qa.add(qa.quantity(decdir), qa.quantity(offset))
        decdir = qa.tos(qa.add(decref, qa.quantity("32arcsec")))
        pcenter = self.__phase_center_string(radir, decdir, dirframe)
        expdir = me.direction(dirframe, radir, decref)
        for myfield in ('1', 'pretend'):
            shift_and_clean(myfield, expdir)
            x = imstat(self.pshift_im + '.image')
            self.assertTrue(
                x['max'][0]/x['rms'][0] < 5,
                msg='Incorrectly found signal in empty field, got S/N of '
                + str(x['max'][0]/x['rms'][0])
            )
            self.__delete_intermediate_products()


class phaseshift_subfunctions_test(unittest.TestCase):

    def test__convert_to_j2000(self):
        from casatasks.private.task_phaseshift import _convert_to_ra_dec_j2000

        phasecenter = 'J2000 19h53m50 40d06m00'
        fra, fdec = _convert_to_ra_dec_j2000(phasecenter)
        places = 6
        self.assertAlmostEqual(fra, -1.074105, places=places)
        self.assertAlmostEqual(fdec, 0.6998770, places=places)

    def test__convert_to_j2000_wrong(self):
        from casatasks.private.task_phaseshift import _convert_to_ra_dec_j2000

        phasecenter = 'BOGUS xxh53m50 40d06m00'
        with self.assertRaisesRegex(RuntimeError, expected_regex="failed"):
            fra, fdec = _convert_to_ra_dec_j2000(phasecenter)


class phaseshift_multi_phasecenter_test(phaseshift_base_checks):
    """ Tests around the use of multi-field phasecenter values (dicts) """

    # Other candidates, could have been:
    # uid___A002_X30a93d_X43e_small.ms: 3 fields, but only 3 scans, and >240MB
    # uid___X02_X3d737_X1_01_small.ms: 3 fields, but only 3 scans
    # twocenteredpointsources.ms: simulated, 2 fields, only 2 scans, (fixvis)
    datadir_multifield = os.path.join('measurementset', 'alma')
    ms_multifield = "uid___A002_X1c6e54_X223-thinned.ms"
    datapath_multifield = ctsys_resolve(os.path.join(datadir_multifield,
                                                     ms_multifield))
    relevant_ddis_multifield = ['0', '1']

    def setUp(self):
        shutil.copytree(self.datapath_multifield, datacopy)
        self.outputvis = "test_vis_multi_field_phasecenter_dict.ms"

    def tearDown(self):
        shutil.rmtree(datacopy)

        if os.path.exists(self.outputvis):
            shutil.rmtree(self.outputvis)

    def get_col_for_field(self, vis, field, ddi, col_name='DATA',):
        """
        Get the visibilities from a data column (specified in param 'col_name'),
        given a field ID and a DDI ID.

        Note: orders the vis values by TIME col, using TaQL, to prevent
        re-ordering seen for example with uid___A002_X1c6e54_X223-thinned.ms
        (It might have needed to reorder by ANTENNA1, ANTENNA2, SCAN_NUMBER,
         etc. but that was not needed)
        """
        tbt = table()

        try:
            tbt.open(vis)

            if field is None:
                col = tbt.getcol(col_name)
            else:
                query_str = ('FIELD_ID in [{}] AND DATA_DESC_ID in [{}]'
                             ' ORDERBY TIME'.format(field, ddi))

                try:
                    query_col = tbt.query(query_str, columns=col_name,
                                          style='python')
                    col = query_col.getcol(col_name)
                finally:
                    query_col.done()
                    if len(col) == 0:
                        raise RuntimeError('Unexpected empty column/query, '
                                           'check test setup')
        finally:
            tbt.done()

        return col

    def check_field_unchanged(self, inputvis, outputvis, field, ddi='0'):
        """ Compare the visibility values between the input and output
        MSs, given one field ID and DDI ID and ensure that the values are
        close (for unchanged / 'passthrough' fields) . """
        vis_in = self.get_col_for_field(datacopy, field=field, ddi=ddi)

        vis_out = self.get_col_for_field(outputvis, field=field, ddi=ddi)

        m1 = np.mean(vis_in)
        m2 = np.mean(vis_out)
        self.assertEqual(m1, m2)
        np.testing.assert_allclose(vis_in, vis_out, rtol=1e-7)

    def test_phasecenter_dict_outofrange(self):
        ''' Check handling of dict with unknown / too many fields '''
        new_center = 'J2000 19h53m50 40d06m00'
        with self.assertRaisesRegex(RuntimeError, "field IDs"):
            result = phaseshift(datacopy, outputvis=self.outputvis,
                                phasecenter={'0': new_center,
                                             '3': new_center})

    def test_phasecenter_dict_simple(self):
        ''' Check multiple field phasecenter(s) given as a dict, one field '''
        new_center = 'J2000 19h53m50 40d06m00'
        phasecenter = {'0': new_center,}
        result = phaseshift(datacopy, outputvis=self.outputvis,
                            phasecenter=phasecenter)

        self.assertEqual(result, None)
        self.check_nrows(self.outputvis, 15344)
        self.check_field_subtable(self.outputvis, datacopy, phasecenter)

        for field in ['1', '2']:
            for ddi in self.relevant_ddis_multifield:
                self.check_field_unchanged(datacopy, self.outputvis, ddi=ddi,
                                           field=field)

    def test_test_nodict(self):
        new_center = 'J2000 19h53m50 40d06m00'
        result = phaseshift(datacopy, outputvis=self.outputvis,
                            phasecenter=new_center)

        self.check_nrows(self.outputvis, 15344)
        self.check_field_subtable(self.outputvis, datacopy, new_center)

    def test_phasecenter_dict_one_out(self):
        ''' Check multiple field phasecenter(s) given as a dict, skip one field '''
        new_centerA = 'J2000 19h53m50 40d06m00'
        new_centerB = 'J2000 22h01m02 40d04m03'
        phasecenter = {'0': new_centerA,
                       '2' : new_centerB}
        result = phaseshift(datacopy, outputvis=self.outputvis,
                            phasecenter=phasecenter)

        self.assertEqual(result, None)
        self.check_nrows(self.outputvis, 15344)
        self.check_field_subtable(self.outputvis, datacopy, phasecenter)

        for ddi in self.relevant_ddis_multifield:
            self.check_field_unchanged(datacopy, self.outputvis, ddi=ddi,
                                       field='1')

    def test_phasecenter_dict_with_field_selection_overlapping(self):
        ''' Check multiple field phasecenter(s) given as a dict,
        skip one field, with selection of all fields in phasecenter dict '''
        new_centerA = 'J2000 19h53m50 40d06m00'
        new_centerB = 'J2000 22h01m02 40d04m03'
        phasecenter = {'0': new_centerA,
                       '2' : new_centerB}
        result = phaseshift(datacopy, outputvis=self.outputvis, field='0,2',
                            phasecenter=phasecenter)

        self.assertEqual(result, None)
        self.check_nrows(self.outputvis, 11536)
        self.check_field_subtable(self.outputvis, datacopy, phasecenter)

    def test_phasecenter_dict_with_field_selection_nonoverlapping(self):
        ''' Check multiple field phasecenter(s) given as a dict,
        skip one field, with selection of fields partially overlapping with
        fields in phasecenter '''
        new_centerA = 'J2000 19h53m50 40d06m00'
        new_centerB = 'J2000 22h01m02 40d04m03'
        phasecenter = {'1': new_centerA,
                       '2' : new_centerB}
        result = phaseshift(datacopy, outputvis=self.outputvis, field='0,2',
                            phasecenter=phasecenter)

        self.assertEqual(result, None)
        self.check_nrows(self.outputvis, 11536)
        self.check_field_subtable(self.outputvis, datacopy, phasecenter)
        for ddi in self.relevant_ddis_multifield:
            self.check_field_unchanged(datacopy, self.outputvis, ddi=ddi,
                                       field='0')


if __name__ == '__main__':
    unittest.main()
