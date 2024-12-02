## CASAtasks

CASAtasks is a self-contained python module that provides the tasks from the [CASA](http://casa.nrao.edu/) project. This package depends on the [casatools](https://open-bitbucket.nrao.edu/projects/CASA/repos/casa6/browse) python module being found in your **PYTHONPATH** at runtime. The casatasks are stateless routines and recipes built on casatools.

## Building CASAtasks

Please refer to the toplevel [readme](../readme.md)

## Available Tasks

| Task Name     | Description                                                                                                                     |
| ------------- | ------------------------------------------------------------------------------------------------------------------------------- |
| accor         | Normalize visibilities based on auto-correlations                                                                               |
| accum         | Accumulate incremental calibration solutions into a calibration table                                                           |
| applycal      | Apply calibrations solutions(s) to data                                                                                         |
| asdmsummary   | Summarized description of an ASDM dataset.                                                                                      |
| bandpass      | Calculates a bandpass calibration solution                                                                                      |
| blcal         | Calculate a baseline-based calibration solution (gain or bandpass)                                                              |
| calstat       | Displays statistical information on a calibration table                                                                         |
| clearcal      | Re-initializes the calibration for a visibility data set                                                                        |
| clearstat     | Clear all autolock locks                                                                                                        |
| concat        | Concatenate several visibility data sets.                                                                                       |
| conjugatevis  | Change the sign of the phases in all visibility columns.                                                                        |
| cvel2         | Regrid an MS or MMS to a new spectral window, channel structure or frame                                                        |
| cvel          | regrid an MS to a new spectral window / channel structure or frame                                                              |
| delmod        | Deletes model representations in the MS                                                                                         |
| exportasdm    | Convert a CASA visibility file (MS) into an ALMA or EVLA Science Data Model                                                     |
| exportfits    | Convert a CASA image to a FITS file                                                                                             |
| exportuvfits  | Convert a CASA visibility data set to a UVFITS file:                                                                            |
| feather       | Combine two images using their Fourier transforms                                                                               |
| fixplanets    | Changes FIELD and SOURCE table entries based on user-provided direction or POINTING table, optionally fixes the UVW coordinates |
| fixvis        | Recalculates (u, v, w) and/or changes Phase Center                                                                              |
| flagcmd       | Flagging task based on batches of flag-commands                                                                                 |
| flagdata      | All-purpose flagging task based on data-selections and flagging modes/algorithms.                                               |
| flagmanager   | Enable list, save, restore, delete and rename flag version files.                                                               |
| fluxscale     | Bootstrap the flux density scale from standard calibrators                                                                      |
| ft            | Insert a source model as a visibility set                                                                                       |
| gaincal       | Determine temporal gains from calibrator observations                                                                           |
| gencal        | Specify Calibration Values of Various Types                                                                                     |
| getantposalma | Retrieve antenna position information from ALMA web service                                                                     |
| getcalmodval  | web service client to retrieve VLA calibrator model data                                                                        |
| getephemtable | Get ephemeris data of a specific ephemeris object from JPL-Horizons system                                                      |
| hanningsmooth | Hanning smooth frequency channel data to remove Gibbs ringing                                                                   |
| imcollapse    | Collapse image along one axis, aggregating pixel values along that axis.                                                        |
| imcontsub     | Estimates and subtracts continuum emission from an image cube                                                                   |
| imdev         | Create an image that can represent the statistical deviations of the input image.                                               |
| imfit         | Fit one or more elliptical Gaussian components on an image region(s)                                                            |
| imhead        | List, get and put image header parameters                                                                                       |
| imhistory     | Retrieve and modify image history                                                                                               |
| immath        | Perform math operations on images                                                                                               |
| immoments     | Compute moments from an image                                                                                                   |
| impbcor       | Construct a primary beam corrected image from an image and a primary beam pattern.                                              |
| importasap    | Convert ASAP Scantable data  into a CASA visibility file (MS)                                                                   |
| importasdm    | Convert an ALMA Science Data Model observation into a CASA visibility file (MS)                                                 |
| importatca    | Import ATCA RPFITS file(s) to a measurement set                                                                                 |
| importfits    | Convert an image FITS file into a CASA image                                                                                    |
| importfitsidi | Convert a FITS-IDI file to a CASA visibility data set                                                                           |
| importgmrt    | Convert a UVFITS file to a CASA visibility data set                                                                             |
| importmiriad  | Convert a Miriad visibility file into a CASA MeasurementSet                                                                     |
| importnro     | Convert NOSTAR data into a CASA visibility file (MS)                                                                            |
| importuvfits  | Convert a UVFITS file to a CASA visibility data set                                                                             |
| importvla     | Import VLA archive file(s) to a measurement set                                                                                 |
| impv          | Construct a position-velocity image by choosing two points in the direction plane.                                              |
| imrebin       | Rebin an image by the specified integer factors                                                                                 |
| imreframe     | Change the frame in which the image reports its spectral values                                                                 |
| imregrid      | regrid an image onto a template image                                                                                           |
| imsmooth      | Smooth an image or portion of an image                                                                                          |
| imstat        | Displays statistical information from an image or image region                                                                  |
| imsubimage    | Create a (sub)image from a region of the image                                                                                  |
| imtrans       | Reorder image axes                                                                                                              |
| imval         | Get the data value(s) and/or mask value in an image.                                                                            |
| initweights   | Initializes weight information in the MS                                                                                        |
| listcal       | List antenna gain solutions                                                                                                     |
| listfits      | List the HDU and typical data rows of a fits file:                                                                              |
| listhistory   | List the processing history of a dataset:                                                                                       |
| listobs       | List the summary of a data set in the logger or in a file                                                                       |
| listpartition | List the summary of a multi-MS data set in the logger or in a file                                                              |
| listsdm       | Lists observation information present in an SDM directory.                                                                      |
| listvis       | List measurement set visibilities.                                                                                              |
| makemask      | Makes and manipulates image masks                                                                                               |
| mstransform   | Split the MS, combine/separate/regrid spws and do channel and time averaging                                                    |
| partition     | Task to produce Multi-MSs using parallelism                                                                                     |
| plotants      | Plot the antenna distribution in the local reference frame:                                                                     |
| plotweather   | Plot elements of the weather table; estimate opacity.                                                                           |
| polcal        | Determine instrumental polarization calibrations                                                                                |
| polfromgain   | Derive linear polarization from gain ratio                                                                                      |
| predictcomp   | Make a component list for a known calibrator                                                                                    |
| rerefant      | Re-apply refant to a caltable                                                                                                   |
| rmfit         | Calculate rotation measure.                                                                                                     |
| rmtables      | Remove tables cleanly, use this instead of rm -rf                                                                               |
| sdbaseline    | Fit/subtract a spectral baseline                                                                                                |
| sdcal         | MS SD calibration task                                                                                                          |
| sdfit         | Fit a spectral line                                                                                                             |
| sdfixscan     | Task for single-dish image processing                                                                                           |
| sdgaincal     | MS SD gain calibration task                                                                                                     |
| sdimaging     | SD task: imaging for total power and spectral data                                                                              |
| sdsmooth      | Smooth spectral data                                                                                                            |
| setjy         | Fills the model column with the visibilities of a calibrator                                                                    |
| simalma       | Simulation task for ALMA                                                                                                        |
| simanalyze    | image and analyze measurement sets created with simobserve                                                                      |
| simobserve    | visibility simulation task                                                                                                      |
| slsearch      | Search a spectral line table.                                                                                                   |
| smoothcal     | Smooth calibration solution(s) derived from one or more sources:                                                                |
| specfit       | Fit 1-dimensional gaussians and/or polynomial models to an image or image region                                                |
| specflux      | Report spectral profile and calculate spectral flux over a user specified region                                                |
| specsmooth    | Smooth an image region in one dimension                                                                                         |
| splattotable  | Convert a downloaded Splatalogue spectral line list to a casa table.                                                            |
| split         | Create a visibility subset from an existing visibility set                                                                      |
| spxfit        | Fit a 1-dimensional model(s) to an image(s) or region for determination of spectral index.                                      |
| statwt        | Compute and set weights based on variance of data.                                                                              |
| tclean        | Radio Interferometric Image Reconstruction                                                                                      |
| uvcontsub     | Continuum fitting and subtraction in the uv plane                                                                               |
| uvmodelfit    | Fit a single component source model to the uv data                                                                              |
| uvsub         | Subtract/add model from/to the corrected visibility data.                                                                       |
| virtualconcat | Concatenate several visibility data sets into a multi-MS                                                                        |
| vishead       | List, summary, get, and put metadata in a measurement set                                                                       |
| visstat       | Displays statistical information from a MeasurementSet, or from a Multi-MS                                                      |
| widebandpbcor | Wideband PB-correction on the output of the MS-MFS algorithm             
