import os
import sys
import numpy as np

from casatasks import casalog
from casatools import ctsys, ms, quanta, calibrater, wvr

def wvrgcal(vis=None, caltable=None, toffset=None, segsource=None,
            sourceflag=None, tie=None, nsol=None, disperse=None, 
	    wvrflag=None, statfield=None, statsource=None, smooth=None,
	    scale=None, spw=None, wvrspw=None,
	    reversespw=None,  cont=None, maxdistm=None,
            minnumants=None, mingoodfrac=None, usefieldtab=None, 
	    refant=None, offsetstable=None, rseed=None):
    """
	Generate a gain table based on Water Vapour Radiometer data.
	Returns a dictionary containing the RMS of the path length variation
	with time towards that antenna (RMS) and the discrepency between the RMS
	path length calculated separately for different WVR channels (Disc.).

	  vis -- Name of input visibility file
											
	              default: none; example: vis='ngc5921.ms'								    

	  caltable -- Name of output gain calibration table
	              default: none; example: caltable='ngc5921.wvr'
								
	  toffset -- Time offset (sec) between interferometric and WVR data

	             default: 0 (ALMA default for cycle 1, for cycle 0 it was -1)

	  segsource -- Do a new coefficient calculation for each source

	             default: True

	  tie -- Prioritise tieing the phase of these sources as well as possible
	         (requires segsource=True)
	             default: [] example: ['3C273,NGC253', 'IC433,3C279']
	
	  sourceflag -- Flag the WVR data for these source(s) as bad and do not produce corrections for it
	               (requires segsource=True)
	               default: [] (none) example: ['3C273']

	  nsol -- Number of solutions for phase correction coefficients during this observation.
	          By default only one set of coefficients is generated for the entire observation.
	          If more sets are requested, then they will be evenly distributed in time throughout'
	          the observation. Values > 1 require segsource=False.
	             default: 1

	  disperse -- Apply correction for dispersion
	             default: False

	  wvrflag -- Regard the WVR data for these antenna(s) as bad and use interpolated values instead
	               default: [] (none) example: ['DV03','DA05','PM02']

	  statfield -- Compute the statistics (Phase RMS, Disc) on this field only
	               default: '' (all)

	  statsource -- Compute the statistics (Phase RMS, Disc) on this source only
	               default: '' (all)

	  smooth -- Smooth the calibration solution on the given timescale
	             default: '' (no smoothing), example: '3s' smooth on a timescale of 3 seconds

	  scale -- Scale the entire phase correction by this factor
	             default: 1. (no scaling)

          spw -- List of the spectral window IDs for which solutions should be saved into the caltable
	             default: [] (all spectral windows), example [17,19,21,23]

          wvrspw -- List of the spectral window IDs from which the WVR data should be taken
	             default: [] (all WVR spectral windows), example [0]

	  reversespw -- Reverse the sign of the correction for the listed SPWs
	                (only needed for early ALMA data before Cycle 0)
	             default: '' (none), example: reversespw='0~2,4'; spectral windows 0,1,2,4

	  cont -- Estimate the continuum (e.g., due to clouds)
                     default: False

          maxdistm -- maximum distance (m) an antenna may have to be considered for being part
	              of the <=3 antenna set for interpolation of a solution for a flagged antenna
		      default: 500

          minnumants -- minimum number of near antennas required for interpolation
	                default: 2

          mingoodfrac -- If the fraction of unflagged data for an antenna is below this value (0. to 1.),
	                 the antenna is flagged.
			 default: 0.8

          usefieldtab -- derive the antenna AZ/EL values from the FIELD rather than the POINTING table
	                 default: False

	  refant -- use the WVR data from this antenna for calculating the dT/dL parameters (can give ranked list)
	              default: '' (use the first good or interpolatable antenna), 
                      examples: 'DA45' - use DA45 
                                ['DA45','DV51'] - use DA45 and if that is not good, use DV51 instead

	  offsetstable -- subtract the temperature offsets in this table from the WVR measurements before
	             using them to calculate the phase corrections
		     default: '' (do not apply any offsets)
		     examples: 'uid___A002_Xabd867_X2277.cloud_offsets' use the given table

          rseed -- set random seed (integer) for the wvrgcal fitting routine to this specific value
                   default: 0 - use internal default value
                   example: 54321

    """

    # make ms tool local 
    myms = ms()
    myqa = quanta()

    ## parameters which are different in format between wvrgcal and wvr.gcal:
    # reverse: only exists in wvr.gcal
    # reversespw: string list in wvrgcal, list wvr.gcal
    # wvrflag: list in wvrgcal, single string in wvr.gcal
    # statfield: single string in wvrgcal, list in wvr.gcal
    # statsource: single string in wvrgcal, list in wvr.gcal
    # refant: list in wvrgcal, list in single string in wvr.gcal
    # rseed: only exists in wvr.gcal
    
    try:
        casalog.origin('wvrgcal')

        if not (type(vis)==str) or not (os.path.exists(vis)):
            raise Exception('Visibility data set not found - please verify the name')

        if (caltable == ""):
            raise Exception("Must provide output calibration table name in parameter caltable.")

        if os.path.exists(caltable):
            raise Exception("Output caltable "+caltable+" already exists - will not overwrite.")

        outdir = os.path.dirname(caltable)
        if outdir == '':
            outdir = '.' 
        if not os.access(outdir, os.W_OK):
            raise Exception("Don't have write permission for output directory "+outdir)
        
        vispar = vis

        smoothpar = 1 # this is for the internal smoothing of wvr.gcal(), which we don't use
        smoothing = -1
        if (type(smooth)==str and smooth!=''):
            smoothing = myqa.convert(myqa.quantity(smooth), 's')['value']
            outputpar = caltable + '_unsmoothed' # as intermediate name before smoothing
        else:
            outputpar = caltable
        
        toffsetpar = toffset

        nsolpar = 1
        if nsol>1:
            if not segsource:
                nsolpar = nsol
            else:
                raise Exception("In order to use nsol>1, segsource must be set to False.")

        segsourcepar = segsource

        sourceflagpar = []
        if segsource and (len(sourceflag)>0):
            sourceflagpar = sourceflag
            for src in sourceflag:
                if not (type(src)==int or type(src)==str) or src=='':
                    raise Exception("List elements of parameter sourceflag must be int or non-empty string.")

        tiepar = []
        if segsource and (len(tie)>0):
            tiepar = tie
            for i in range(0,len(tie)):
                src = tie[i]
                if not (type(src)==str) or src=='':
                    raise Exception("List elements of parameter tie must be non-emptystrings.")
                                                
        spwpar = spw
        if (len(spw)>0):
            for myspw in spw:
                if not (type(myspw)==int) or myspw<0:
                    raise Exception("List elements of parameter spw must be int and >=0.")

        wvrspwpar = wvrspw
        if (len(wvrspw)>0):
            for myspw in wvrspw:
                if not (type(myspw)==int) or myspw<0:
                    raise Exception("List elements of parameter wvrspw must be int and >=0.")
                        
        reversespwpar = []
        if not (reversespw==''):
            spws = myms.msseltoindex(vis=vis,spw=reversespw)['spw']
            for id in spws:
                reversespwpar.append(id)

        dispersepar = disperse                
        if disperse:
            dispdirpath = os.getenv('WVRGCAL_DISPDIR', '')
            if not os.path.exists(dispdirpath+'/libair-ddefault.csv'):
                path1 = dispdirpath
                dispdirpath = ctsys.resolve("alma/wvrgcal")
                if not os.path.exists(dispdirpath+'/libair-ddefault.csv'):
                    raise Exception("Dispersion table libair-ddefault.csv not found in path "\
                        +"given by WVRGCAL_DISPDIR nor in \""+dispdirpath+"\"")

                os.putenv('WVRGCAL_DISPDIR', dispdirpath)
                                
            casalog.post('Using dispersion table '+dispdirpath+'/libair-ddefault.csv')

        contpar = cont
        if cont and segsource:
            raise Exception("cont and segsource are not permitted to be True at the same time.")

        usefieldtabpar = usefieldtab
        
        offsetspar = offsetstable

        wvrflagpar = ""
        if (len(wvrflag)>0):
            for ant in wvrflag:
                if not (type(ant)==int or type(ant)==str):
                    raise Exception("List elements of parameter wvrflag must be int or string.")
                if (ant != ''):
                    if len(wvrflagpar)>0:
                        wvrflagpar += ","
                    wvrflagpar += str(ant)

        refantpar = ""
        if (type(refant)!=list):
            refant = [refant]
        if (len(refant)>0):
            for ant in refant:
                if not (type(ant)==int or type(ant)==str):
                    raise Exception("Parameter refant must be int or string or a list of them.")
                if (ant != ''):
                    if len(refantpar)>0:
                        refantpar += ","
                    refantpar += str(ant)

        statfieldpar = []
        if not (statfield==None or statfield=="") and type(statfield)==str:
            statfieldpar = [statfield]

        statsourcepar = []
        if not (statsource==None or statsource=="") and type(statsource)==str:
            statsourcepar = [statsource]

        scalepar = scale
        
        maxdistmpar = maxdistm
                
        minnumantspar = minnumants

        mingoodfracpar = mingoodfrac

        rseedpar = 0
        if type(rseed)==int and rseed>=0:
            rseedpar = rseed
        elif not (rseed==None or rseed==""):
            raise Exception("Parameter rseed must be an integer >= 0 (the value 0 will use the internal default seed).")
            
        casalog.post('Running wvr.gcal ...')


        templogfile = 'wvrgcal_tmp_'+str(np.random.randint(1E6,1E8))
        if not os.access(".", os.W_OK):
            import tempfile
            templogfile = tempfile.gettempdir()+"/"+templogfile

        os.system('rm -rf '+templogfile)

        mywvr = wvr()

        rval = mywvr.gcal(vis=vispar,
                          output=outputpar,
                          toffset=toffsetpar,
                          nsol=nsolpar,
                          segsource=segsourcepar,
                          reverse=False, # only reverse those SPWs in reversespwpar
                          reversespw=reversespwpar,
                          disperse=dispersepar,
                          cont=contpar,
                          wvrflag=wvrflagpar,
                          sourceflag=sourceflagpar,
                          statfield=statfieldpar,
                          statsource=statsourcepar,
                          tie=tiepar,
                          smooth=smoothpar,
                          scale=scalepar,
                          maxdistm=maxdistmpar,
                          minnumants=minnumantspar,
                          mingoodfrac=mingoodfracpar,
                          usefieldtab=usefieldtabpar,
                          spw=spwpar,
                          wvrspw=wvrspwpar,
                          refant=refantpar,
                          offsets=offsetspar,
                          rseed=rseedpar,
                          logfile=templogfile)

        loglines = []
        with open(templogfile) as f:
            loglines = f.readlines()
        for i in range(len(loglines)):
            loglines[i] = loglines[i].expandtabs()
        casalog.post(''.join(loglines))
        
        # prepare variables for parsing log lines to extract info table
        hfound = False
        hend = False
        namel = []
        wvrl = []
        flagl = []
        rmsl = []
        discl = []
        flagsd = {}
        parsingok = True
                
        for ll in loglines:
            if hfound:
                if "Expected performance" in ll:
                    hend = True
                elif not hend:
                    vals = ll.split()
                    wvrv = False
                    flagv = False
                    rmsv = 0.
                    discv = 0.
                    if(len(vals)!=6):
                        casalog.post('Error parsing wvrgcal info table.line: '+ll,'WARN')
                        parsingok=False
                    else:
                        if vals[2]=='Yes':
                            wvrv=True
                        else:
                            wvrv=False
                        if vals[3]=='Yes':
                            flagv=True
                        else:
                            flagv=False
                    try:
                        rmsv = float(vals[4])
                    except:
                        casalog.post('Error parsing RMS value in info table line: '+ll,'WARN')
                        rmsv = -1.
                        parsingok=False
                    try:
                        discv = float(vals[5])
                    except:
                        casalog.post('Error parsing Disc. value in info table line: '+ll,'WARN')
                        discv = -1.
                        parsingok=False

                    namel.append(vals[1])
                    wvrl.append(wvrv)
                    flagl.append(flagv)
                    rmsl.append(rmsv)
                    discl.append(discv)

                                        
            elif (rval==0) and (not hend) and ("Disc (um)" in ll):
                hfound = True
            elif 'WVR data points for antenna' in ll: # take note of antennas flagged because of too few good WVR data
                token = ll.split('for antenna ')[1].split()
                antennaID = int(token[0])
                if 'All WVR' in ll:
                    flagsd[antennaID] = 0.
                else:
                    unflagged = int(token[2])
                    total = float(token[5]) # ends in a period
                    if total>0:
                        flagsd[antennaID] = unflagged / total
                    else:
                        casalog.post('Error: zero datapoints reported for antenna id '+str(antennaID)+' in info table line: '+ll,'WARN')
                        parsingok=False
                        
        # end for ll

        # create list of flagging fractions for each antenna
        unflagfracl = list(np.ones(len(namel)))
        for myid in flagsd.keys():
            if myid >= len(unflagfracl):
                casalog.post('Error: flagged antenna id '+str(myid)+' > max known antenna id '+str(len(unflagfracl)-1) ,'WARN')
                parsingok=False
            else:
                unflagfracl[myid] = flagsd[myid]
                    

        os.system('rm -rf '+templogfile)
        
        taskrval = { 'Name': namel,
                     'WVR': wvrl,
                     'Flag': flagl,
                     'Frac_unflagged': unflagfracl,
                     'RMS_um': rmsl,
                     'Disc_um': discl,
                     'rval': rval,
                     'success': False}
        
        for k in range(len(namel)):
            if(flagl[k] and rmsl[k]==0. and discl[k]==0.):
                casalog.post('Solution for flagged antenna '+namel[k]
                             +' could not be interpolated due to insufficient number of near antennas. Was set to unity.',
                             'WARN')

        if (rval==0) and parsingok:
            taskrval['success'] = True
                

        if(rval == 0):
            if (smoothing>0):
                mycb = calibrater()
                mycb.open(filename=vis, compress=False, addcorr=False, addmodel=False)
                mycb.smooth(tablein=caltable+'_unsmoothed', tableout=caltable,
                            smoothtype='mean', smoothtime=smoothing)
                mycb.close()
            return taskrval
        else:
            if(rval == 255):
                casalog.post('wvr.gcal terminated with exit status '+str(rval),'SEVERE')
                return taskrval
            elif(rval == 134 or rval==1):
                casalog.post('wvr.gcal terminated with exit status '+str(rval),'WARN')
                casalog.post("No useful input data.",'SEVERE')
                return taskrval
            else:
                casalog.post('wvrgcal terminated with exit status '+str(rval),'WARN')
                return taskrval
        
    except Exception as instance:
        print('*** Error *** ', instance)
        raise Exception from instance
