from __future__ import absolute_import

import os
import numpy as np

# Assume CASA 6, CASA 5 is no longer built
from casatools import image
ia = image()


class ImagingDict():
    """
    Keeps track of tclean + deconvolve return dictionaries
    """
    
    def __init__(self) -> None:
        self._returndict = {}
        self._summaryminor = {}
        
    def __str__(self) -> str:
        """
        Pretty print the return dictionary.
        """
        
        retstr = ''
        for key, val in self._returndict.items():
            if key == 'summaryminor':
                retstr += 'summaryminor:\n'
                for field in val.keys():
                    for freq in val[field].keys():
                        for stokes in val[field][freq].keys():
                            for key2, val2 in val[field][freq][stokes].items():
                                retstr += f'\t Field: {field: 3d}, Chan: {freq: 3d}, Stokes: {stokes :2d}, {key2}, {val2}\n'
            else:
                retstr += f'{key}: {val}\n'
        return retstr
        
    @property
    def returndict(self) -> dict:
        return self._returndict
    
    @returndict.setter
    def returndict(self, inpdict: dict) -> None:
        self._returndict = inpdict
        
    @returndict.deleter
    def returndict(self) -> None:
        del self._returndict


    def initialize_dict(self) -> None:
        """
        Initialize the return dictionary with dummy values. This sets up the
        dictionary structure and keys to expect, and will be filled in with
        real values when running tclean or deconvolve.
        """
        
        # Initialize the values that don't need to inspect any images
        self._returndict['cleanstate'] = ''
        self._returndict['cyclefactor'] = 1.0
        self._returndict['cycleiterdone'] = 0
        self._returndict['cycleniter'] = 0
        self._returndict['cyclethreshold'] = 0
        self._returndict['interactiveiterdone'] = 0
        self._returndict['interactivemode'] = False
        self._returndict['interactiveniter'] = 0
        self._returndict['interactivethreshold'] = 0

        self._returndict['iterdone'] = 0
        self._returndict['loopgain'] = 0
        self._returndict['maxpsffraction'] = 0
        self._returndict['maxpsfsidelobe'] = 0
        self._returndict['minpsffraction'] = 0

        self._returndict['niter'] = 0
        self._returndict['nmajordone'] = 0
        self._returndict['nsigma'] = 0.0
        self._returndict['stopcode'] = 0

        self._returndict['summarymajor'] = np.array([])
        # Summary minor is nested as {field{freq{stokes}}
        self._returndict['summaryminor'] = {}
        self._returndict['summaryminor'][0] = {}
        self._returndict['summaryminor'][0][0] = {}
        self._returndict['summaryminor'][0][0][0] = self.initialize_summary_minor()

        self._returndict['threshold'] = 0.0
        self._returndict['stopDescription'] = 'Zero iterations performed'
                                      

    def initialize_summary_minor(self) -> dict:
        """
        Initialize the summary minor dictionary with dummy values. This sets up the
        dictionary structure and keys to expect, and will be filled in with
        real values when running tclean or deconvolve.
        """

        self._summaryminor['iterDone'] = []
        self._summaryminor['peakRes'] = []
        self._summaryminor['modelFlux'] = []
        self._summaryminor['cycleThresh'] = []
        self._summaryminor['cycleStartIters'] = []
        self._summaryminor['startIterDone'] = []
        self._summaryminor['startPeakRes'] = []
        self._summaryminor['startModelFlux'] = []
        self._summaryminor['startPeakResNM'] = []
        self._summaryminor['peakResNM'] = []
        self._summaryminor['masksum'] = []
        self._summaryminor['mpiServer'] = []
        self._summaryminor['stopCode'] = []
                                      
        return self._summaryminor
                                      


    def get_key(self, key:str, field:int=0, chan:int=0, stokes:int=0) -> list:
        """
        Return the list of values for the specified key. If requesting a key
        within "summaryminor", also specify the field, channel and stokes plane.
                                      
        Inputs:
        key         The key to return. str
        field       The field to return. int
        chan        The channel to return. int
        stokes      The stokes plane to return. int

        Returns:
        The list of values for the specified key. list
        """
                                      
        summaryminor_keys = ['iterDone', 'peakRes', 'modelFlux', 'cycleThresh',
                             'cycleStartIters', 'startIterDone',
                             'startPeakRes', 'startModelFlux',
                             'startPeakResNM', 'peakResNM', 'masksum',
                             'mpiServer', 'stopCode']
        
        try:
            if key in summaryminor_keys:
                return self._returndict['summaryminor'][field][chan][stokes][key]
            else:
                return self._returndict[key]
        except KeyError:
            print('WARNING : Key not found in return dictionary.')
            return []
                                      

    def append_dict(self, inpdict:dict) -> None:
        """
        Append the input dictionary to the return dictionary. This is used to
        merge dictionaries from multiple tclean/deconvolve calls.
        
        The input dictionary must be a fully formed return dictionary from
        either `tclean` or `deconvolve`.
                                      
        The returndict attribute is modified in place.
        
        Inputs:
        inpdict     The input dictionary to append. dict
                                      
        Returns:
        None
        """
                                      
        # Only append these keys, replace the rest
        append_keys = ['summarymajor', 'summaryminor']
        # Accumulate these keys, incrementing with every append
        accum_keys = ['cycleniter', 'cycleiterdone', 'interactiveiterdone', 'interactiveniter', 'iterdone', 'nmajordone', 'niter']
        
        try:
            for key, val in inpdict.items():
                if key in append_keys:
                    if key == 'summarymajor':
                        self.returndict[key] = np.append(self.returndict[key], val)
                    elif key == 'summaryminor':
                        # Iterate through the dict and append individual values
                        for field in val.keys():
                            for freq in val[field].keys():
                                for stokes in val[field][freq].keys():
                                    for key2, val2 in val[field][freq][stokes].items():
                                        self._returndict[key][field][freq][stokes][key2] = np.append(self._returndict[key][field][freq][stokes][key2], val2)
                elif key in accum_keys:
                    self._returndict[key] += val
                else:
                    self._returndict[key] = val
        except KeyError:
            print('Input dictionary does not have expected keys.')
            raise


class ReturnDictionary():
    """
    Class that constructs the tclean return dictionary - implemented
    specifically for the niter=0 case, to reproduce the dictionary generated in
    C++ without modifying the C++ code.

    Attributes:
    -----------
    residname       : Name of the input residual image (required)
    modelname       : Name of the input model image (optional)
    summaryminor    : The dictionary containing summaryminor
    retrec          : The full return dictionary

    Methods:
    --------

    imageDimensions(residname)
        Use the image() tool to query for the number of Stokes planes and
        frequency channels in the input image.

    fillSummaryMinor(residname, modelname, channo, stokes, stokes_axis, freq_axis, fullsummary):
        Given the input image name, and the corresponding field, channel number, and
        Stokes plane, extract the relevant information from the image to generate a
        summaryminor dict as defined

    constructSummaryMinor(self, paramList):
        Constructs and populates a nested dictionary containing the summaryMinor()
        information required by tclean.
        This is a duplicate of the existing summaryMinor() class in C++, but is
        being re-written here to avoid issues with modifying iterControl.

    constructResidualDict(paramList):
        Construct the residual dictionary given the input image name and associated parameters. This residual dictionary is meant to
        duplicate that generated by imager.getSummary(fullsummary) for the special
        case of niter = 0 to avoid initializing the deconvolver.
    """

    def __init__(self):
        """
        Initialize necessary variables for the return dictionary.
        """
        self.residname = ''
        self.modelname = ''
        self.summaryminor = dict()
        self.retrec = dict()


    def imageDimensions(self):
        """
        Given the input image, uses the ia tool to query for the number of Stokes
        planes and frequency channels.

        Inputs:
        None

        Returns:
        nstokes         Number of Stokes planes in the image
        nfreq           Number of frequency channels in the image
        """
        ia = image()

        # At this point we've already checked residname exists
        ia.open(self.residname)
        csys = ia.coordsys()

        shape = ia.shape()

        #  Figure out which axis is which
        stokes_axis = csys.findaxisbyname('Stokes')
        freq_axis = csys.findaxisbyname('Frequency')

        nstokes = shape[stokes_axis]
        nfreq = shape[freq_axis]

        ia.close()

        return nstokes, nfreq, stokes_axis, freq_axis


    def fillSummaryMinor(self, channo, stokes, stokes_axis, freq_axis, fullsummary):
        """
        Given the input image name, and the corresponding field, channel number, and
        Stokes plane, extract the relevant information from the image to generate a
        summaryminor dict as defined
        https://casadocs.readthedocs.io/en/stable/notebooks/synthesis_imaging.html#Minor-Cycle-Summary-Dictionary

        Inputs:
        channo          Channel number to query in the image, int
        stokes          Stokes plane to query in the image, int
        stokes_axis     The axis to index for Stokes, int
        freq_axis       The axis to index for frequency, int
        fullsummary     Construct a full summary, or only a subset, bool

        Returns:
        summaryparams     Dict containing the necessary (key:value) pairs
        """

        ia = image()

        if not os.path.exists(self.residname):
            raise FileNotFoundError(f'Residual image {self.residname} does not exist.')

        ia.open(self.residname)
        shape = ia.shape()

        if stokes_axis == 2 and freq_axis == 3:
            blc = [0, 0, stokes, channo]
            trc = [shape[0], shape[1], stokes, channo]
        elif stokes_axis == 3 and freq_axis == 2:
            blc = [0, 0, channo, stokes]
            trc = [shape[0], shape[1], channo, stokes]

        data = ia.getchunk(blc, trc, dropdeg=True)
        mask = ia.getchunk(blc, trc, dropdeg=True, getmask=True)
        ia.close()

        # If model image exists, calc model flux, else set to 0
        model_sum = 0
        if os.path.exists(self.modelname):
            ia.open(self.modelname)
            model_data = ia.getchunk(blc, trc, dropdeg=True)
            model_sum = np.sum(model_data)
            ia.close()

        peak_resid = np.amax(data*mask)
        if fullsummary:
            peak_resid_NM = np.amax(data)
            mask_sum = np.sum(mask)

        summaryparams = dict()
        # This entire function is only invoked in the special case of niter=0
        summaryparams['iterDone'] = [0.0,]
        summaryparams['peakRes'] = [peak_resid,]
        # model flux has to be zero because no iterations were performed
        summaryparams['modelFlux'] = [model_sum,]
        # No threshold because no deconvolution done
        summaryparams['cycleThresh'] = [0.0,]

        if fullsummary:
            summaryparams['cycleStartIters'] = [0.0,]
            summaryparams['startIterDone'] = [0.0,]
            summaryparams['startPeakRes'] = [peak_resid,]
            summaryparams['startModelFlux'] = [model_sum,]
            summaryparams['startPeakResNM'] = [peak_resid_NM,]
            summaryparams['peakResNM'] = [peak_resid_NM,]
            summaryparams['masksum'] = [mask_sum,]
            summaryparams['mpiServer'] = [0.0,]
            summaryparams['stopCode'] = [3,]


        return summaryparams


    def constructSummaryMinor(self, paramList):
        """
        Constructs and populates a nested dictionary containing the summaryMinor()
        information required by tclean.

        This is a duplicate of the existing summaryMinor() class in C++, but is
        being re-written here to avoid issues with modifying iterControl.

        Inputs:
        paramList        Object that contains all imaging parameters

        Returns:
        summaryMinor    Nested dictionary with the summaryMinor keys, dict
        """

        impars = paramList.allimpars
        decpars = paramList.alldecpars

        # Each field is stored as a different key in impars
        nfields = len(impars.keys())

        for ff in range(nfields):
            self.residname=impars[str(ff)]['imagename']+'.residual.tt0' if(os.path.exists(impars[str(ff)]['imagename']+'.residual.tt0')) else impars[str(ff)]['imagename']+'.residual'
            self.modelname=impars[str(ff)]['imagename']+'.model.tt0' if(os.path.exists(impars[str(ff)]['imagename']+'.model.tt0')) else impars[str(ff)]['imagename']+'.model'

            fullsummary = decpars[str(ff)]['fullsummary']
            nstokes, nfreq, stokes_axis, freq_axis = self.imageDimensions()

            self.summaryminor[ff] = dict()

            for cc in range(nfreq):
                self.summaryminor[ff][cc] = dict()
                for ss in range(nstokes):
                    self.summaryminor[ff][cc][ss] = self.fillSummaryMinor(cc, ss, stokes_axis, freq_axis, fullsummary)

        return self.summaryminor



    def constructResidualDict(self, paramList):
        """
        Construct the residual dictionary given the input image name and associated parameters. This residual dictionary is meant to
        duplicate that generated by imager.getSummary(fullsummary) for the special
        case of niter = 0 to avoid initializing the deconvolver.

        Inputs:
        paramList   The tclean inputs, defaults where not specified. dict

        Returns:
        retrec      The return dictionary with imaging statistics. dict
        """

        imagename = paramList.allimpars['0']['imagename']

        self.residname = imagename+'.residual.tt0' if(os.path.exists(imagename +'.residual.tt0')) else imagename +'.residual'
        do_summary_minor = True
        # Only fill summaryminor if the residual image exists
        if not os.path.exists(self.residname):
            do_summary_minor = False
            #raise FileNotFoundError(f'{residname} does not exist on disk. Cannot construct tclean return dictionary.')

        # Initialize the values that don't need to inspect any images
        self.retrec['cleanstate'] = 'running'
        self.retrec['cyclefactor'] = paramList.getAllPars()['cyclefactor']
        self.retrec['cycleiterdone'] = 0
        self.retrec['cycleniter'] = 0
        self.retrec['cyclethreshold'] = 0
        self.retrec['interactiveiterdone'] = 0
        self.retrec['interactivemode'] = paramList.alldecpars['0']['interactive']
        self.retrec['interactiveniter'] = 0
        self.retrec['interactivethreshold'] = 0

        self.retrec['iterdone'] = 0
        self.retrec['loopgain'] = 0
        self.retrec['maxpsffraction'] = 0
        self.retrec['maxpsfsidelobe'] = 0
        self.retrec['minpsffraction'] = 0

        self.retrec['niter'] = 0
        self.retrec['nmajordone'] = 1
        self.retrec['nsigma'] = 0.0
        # stopcode 3 --> Zero iterations performed
        self.retrec['stopcode'] = 3

        self.retrec['summarymajor'] = np.array([0,])

        if do_summary_minor:
            self.retrec['summaryminor'] = self.constructSummaryMinor(paramList)
        else:
            self.retrec['summaryminor'] = {}

        self.retrec['threshold'] = 0.0
        self.retrec['stopDescription'] = 'Zero iterations performed'


        return self.retrec
