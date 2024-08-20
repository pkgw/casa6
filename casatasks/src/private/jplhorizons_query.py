'''
 main module for retreiving JPL-Horizons ephemeris data and
 convert it to CASA readable format (i.e. CASA table )
'''
import os
import time
from math import sqrt, sin, cos, log
import numpy as np

from casatools import table, quanta, measures
from casatasks import casalog

_tb = table()
_qa = quanta()
_me = measures()

# debug
_debug = False 

def gethorizonsephem(objectname, starttime, stoptime, incr, outtable, asis=False, rawdatafile=''):
    """
    Main driver function for ephemeris data query from JPL-Horizons
    
    arguments:
        objectname: ephemeris object name (case insensitive). Try to convert 
                    a common name to ID if asis = False.
        starttime: start time of the ephemeris data (expects YYYY/MM/DD/HH:MM 
                   format
        stoptime: stop (end) time of the ephemeris data in the same format as
                  starttime
        incr: time increment (interval) of the ephemeris data
        outtable: output CASA table name
        asis: a toggle for additinal check to resolve object name 
        rawdatafile: raw query result file name (optional)
    """

    # commented ones are not currently supported in setjy
    # Juno's table is exist in the data repo but not supported in setjy
    asteroids = {'ceres': '1;',
                 'pallas': '2;',
                 'juno': '3;',  # large crater and temperature varies
                 'vesta': '4;',
                 'lutetia': '21;'
                 }
    planets_and_moons = {'sun': '10',
                         'mercury': '199',
                         'venus': '299',
                         'moon': '301',
                         'mars': '499',
                         'jupiter': '599',
                         'io': '501',
                         'europa': '502',
                         'ganymede': '503',
                         'callisto': '504',
                         'saturn': '699',
                         'titan': '606',
                         'uranus': '799',
                         'neptune': '899',
                         'pluto': '999'}
    known_objects = planets_and_moons
    known_objects.update(asteroids)
    # default quantities (required by setjy) for CASA
    quantities = '1,14,15,17,19,20,24'
    ang_format = 'DEG'

    # check objectname. If it matches with the known object list. If asis is specified
    # it will skip the check.
    # For small bodies such as comets, objectname must include 'DES=' following designated
    # object name or id (SPK ID). Based on the JPL-Horizons documentaiton, for IAU ID, 'DES='
    # can be omitted. https://ssd.jpl.nasa.gov/horizons/app.html#command
    start_time = None
    stop_time = None
    step_size = None
    if not asis:
        if not objectname.lower() in known_objects:
            raise ValueError(
                f"{objectname} is not in the known object list for CASA. To skip this check set asis=True")
        else:
            target = known_objects[objectname.lower()]
    else:
        target = objectname
    try:
        if starttime.startswith('JD') or starttime.startswith('MJD'):
            start_time = starttime
            stop_time = stoptime
        else:
            start_time = _qa.time(starttime, form='ymd')
            stop_time = _qa.time(stoptime, form='ymd')
    except ValueError as e:
        casalog.post(e)

    try:
        step_size = incr.replace(' ', '')
    except ValueError as e:
        casalog.post(e)
   
    if _debug:
        print("target=",target)
        print(f"start_time={start_time}, stop_time={stop_time}, step_size={step_size}")
        print(f"quantities={quantities}, ang_format={ang_format}")
    ephemdata = queryhorizons(target, start_time, stop_time, step_size, quantities, ang_format, rawdatafile)
    # return ephemdata
    if ephemdata and 'result' in ephemdata:
        casalog.post('converting ephemeris data to a CASA table')
        tocasatb(ephemdata, outtable)
    if not os.path.exists(outtable):
        casalog.post('Failed to produce the CASA ephemeris table', 'WARN') 
        if rawdatafile=='':
            casalog.post('You may want to re-run the task with rawdatafile'+\
                ' parameter specified and check the content of the output raw data file for the further information on the error.','WARN') 
        else:
            casalog.post(f'Please check the content of {rawdatafile} for the further information of the error','WARN')
               

def queryhorizons(target, starttime, stoptime, stepsize, quantities, ang_format, rawdatafile=''):
    """
    Submit a query to the JPL-Horizons API
   
    arguments:
        target: a target solar system object name/id (str)
        starttime: ephemeris start time (e.g. '2021/12/01/00:00')
        stoptime: ephemeris stop time
        stepsize: ephemeris data time increment
        quantities: output data quantities given as a list of integers which represent the specific 
                    data as defined in https://ssd,jpl.nasa.gov/horizons/manual.html
        ang_format: RA, Dec output format ('DEG' or 'HMS') 
        rawdatafile: a file name to save raw query results (Optional)

    """
    import ast
    import certifi
    import ssl
    from urllib.request import Request, urlopen
    from urllib.parse import urlencode
    from urllib.error import URLError
    import socket

    context = ssl.create_default_context(cafile=certifi.where())
    urlbase = 'https://ssd.jpl.nasa.gov/api/horizons.api'
    data = None

    values = {'format': 'json',
              'EPHEM_TYPE': 'OBSERVER',
              'OBJ_DATA': 'YES',
              'COMMAND': f"'{target}'",
              'START_TIME': starttime,
              'STOP_TIME': stoptime,
              'STEP_SIZE': stepsize,
              'CENTER': '500@399',
              'QUANTITIES': f"'{quantities}'",
              'ANG_FORMAT': ang_format}

    pardata = urlencode(values, doseq=True, encoding='utf-8')
    params = pardata.encode('ascii')
    # for debugging
    #print("params=", params)

    req = Request(urlbase, params)
    datastr = None
    try:
        with urlopen(req, context=context, timeout=10) as response:
            datastr = response.read().decode()
    except URLError as e:
        if hasattr(e, 'reason'):
            msg = f'URLError: Failed to reach URL={urlbase} (reason: {e.reason})'
            casalog.post(msg, 'WARN')
        if hasattr(e, 'code'):
            msg = f'URLError: Couldn\'t fullfill the request to {urlbase} (code: {e.code})'
            casalog.post(msg, 'WARN')
    except socket.timeout as e:
        msg = f'Failed to reach URL={urlbase}. Socket timeout {e}'
        casalog.post(msg, 'WARN')

    status = response.getcode()
    if datastr:
        try:
            data = ast.literal_eval(datastr)
        except RuntimeError as e:
            casalog.post(e, 'SEVERE')
        if status == 200:
           #
           # If the ephemeris data file was generated, write it to the output file:
            if rawdatafile != "":
                if os.path.exists(rawdatafile):
                    os.remove(rawdatafile)
                if "result" in data:
                    with open(rawdatafile, "w") as f:
                        f.write(data["result"])
                # Otherwise, the ephemeris file was not generated so output an error
                else:
                    casalog.post("ERROR: No data found. Ephemeris data file not generated", 'WARN')
        else:
            raise RuntimeError(f'Could not retrieve the data. Error code:{status}:{response.msg}')
    else:
        data = None
    return data


def tocasatb(indata, outtable):
    """
    convert a JPL-Horizons query results to a CASA table
    indata: either query data ('results') or file name that contains the result data
    outtable: output CASA table name
    """
    # input data  columns 
    # Date__(UT)__HR:MN     R.A.___(ICRF)___DEC  ObsSub-LON ObsSub-LAT  SunSub-LON SunSub-LAT ¥
    # NP.ang   NP.dist                r        rdot             delta      deldot     S-T-O
    import re

    # output CASA table columns     
    # These are required columns for SetJy
    cols = {
        'MJD': {'header': 'Date__\(UT\)',
                'comment': 'date in MJD',
                'unit': 'd'},
        'RA': {'header': 'R.A.',
               'comment': 'astrometric Right Ascension (ICRF)',
               'unit': 'deg'},
        'DEC': {'header': 'DEC',
                'comment': 'astrometric Declination (ICRF)',
                'unit': 'deg'},
        'Rho': {'header': 'delta',
                'comment': 'geocentric distance',
                'unit': 'AU'},
        'RadVel': {'header': 'deldot',
                   'comment': 'geocentric distance rate',
                   'unit': 'AU/d'},
        'NP_ang': {'header': 'NP.ang',
                   'comment': 'North-Pole pos. angle',
                   'unit': 'deg'},
        'NP_dist': {'header': 'NP.dist',
                    'comment': 'North-Pole ang. distance',
                    'unit': 'arcsec'},
        'DiskLong': {'header': 'ObsSub-LON',
                     'comment': 'Sub-observer longitude',
                     'unit': 'deg'},
        'DiskLat': {'header': r'ObsSub-LAT',
                    'comment': 'Sub-observer latitude',
                    'unit': 'deg'},
        'Sl_lon': {'header': 'SunSub-LON',
                   'comment': 'Sub-Solar longitude',
                   'unit': 'deg'},
        'Sl_lat': {'header': r'SunSub-LAT',
                   'comment': 'Sub-Solar longitude',
                   'unit': 'deg'},
        'r': {'header': 'r',
              'comment': 'heliocentric distance',
              'unit': 'AU'},
        'rdot': {'header': 'rdot',
                 'comment': 'heliocentric distance rate',
                 'unit': 'km/s'},
        'phang': {'header': 'S-T-O',
                  'comment': 'phase angle',
                  'unit': 'deg'}
    }

    # do in a two-step
    # read the original query result dictionary containing ephemeris data and save the data part
    # to a temporary text file. Then re-read the temporary file to a casa table
    # with the approriate data format that setjy and measure expect.
    tempfname = 'temp_ephem_'+str(os.getpid())+'.dat'
    tempconvfname = 'temp_ephem_conv_'+str(os.getpid())+'.dat'
    try:
        exceedthelinelimit = None
        ambiguousname = None
        queryerrmsg = ''
        ephemdata = None
        # Scan the original data
        if isinstance(indata, dict) and 'result' in indata:
            #print("ephem data dict")
            ephemdata = indata['result']
        elif isinstance(indata, str):
            if os.path.exists(indata):
                with open(indata, 'r') as infile:
                    #ephemdata = infile.readlines()
                    ephemdata = infile.read()
                    # the relevant information in the main header section is extracted as
            else:
             raise IOError(f'{indata} does not exist')
        else:
             raise IOError(f'{indata} should be a dictionary contains ephemeris data under'+
                            ' the keyword, "results" or the text file name contains the raw JPL-Horizons'+
                            ' query results')
        # it is read line by line. The ephemeris data is stored in the separate text file
        # to be further processed in subsequent steps.
        with open(tempfname, 'w') as outfile:
            # Some initialization
            headerdict = {}
            # jplhorizonsdataIdFound = False
            datalines = 0
            readthedata = False
            startmjd = None
            endmjd = None
            incolnames = None
            # multiple entries (in different units) may exit for orb. per.
            foundorbper = False
            ###
            lcnt = 0
            # for lnum, line in enumerate(infile):
            for lnum, line in enumerate(ephemdata.split('\n')):
                # JPL-Horizons data should contain this line at the beginning
                if re.search(r'JPL/HORIZONS', line):
                    # jplhorizondataIdFound = True
                    casalog.post("Looks like JPL-Horizons data","INFO2")
                elif re.search(r'^\s*Ephemeris\s+', line):  # date the data file was retrieved and  created
                    #m = re.search(r'(API_USER\S+\s+(\S+)\s+([0-9]+)\s+(\S+)\s+(\S+)')
                    (_, _, _, wday, mon, day, tm, year, _) = re.split(' +', line, 8)
                    # date info in 3-7th items

                    try:
                        float(mon)
                        nmon = mon
                    except:
                        tonummon=time.strptime(mon,"%b")
                        nmon = f"{tonummon.tm_mon:02d}"
                    day2 = f"{int(day):02d}"
                    headerdict['VS_CREATE'] = year + '/' + nmon + '/' + day2 + '/' + tm[0:5]
                    # VS_DATE - use the current time to indicate the time CASA table is created
                    headerdict['VS_DATE'] = time.strftime('%Y/%m/%d/%H:%M',time.gmtime() )
                    headerdict['VS_TYPE'] = 'Table of comet/planetary positions'
                    # VERSION stored in the output table may be incremented in the future.
                    # For now, it is fixed, but it may be incremented from 0003 to 0004 to indiate
                    # this new code is used to convert the jpl horizons data to a table.
                    # ver 0004.0001 - new keyword: ephemeris_source, posrefsys 'ICRF/J2000' -> "ICRS' 
                    headerdict['VS_VERSION'] = '0004.0001'
                    # target object name
                elif re.match(r'^[>\s]*Target body name', line):
                    m = re.match(r'^[>\s]*Target body name:\s+(\S+)\s+(\S*)\s+\{(\w+)\:\s*(\S+)\}', line)
                    if m:
                        headerdict['NAME'] = m[1]
                        if len(m.groups())==4 and m[3]=='source':
                            headerdict['ephemeris_source'] = m[4]
                # start time (of the requested time range)
                elif re.search(r'Start time', line):
                    m = re.match(r'^[>\s]*Start time\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\w+)', line)
                    if m:
                        startmjd = _qa.totime(m[1] + '/' + m[2])
   #--This info will not be used but left here since it might be useful fo debugging.
   #             # end time (of the requested time range)
                elif re.search(r'Stop  time', line):
                    m = re.match(r'^[>\s]*Stop  time\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\w+)', line)
                    if m:
                        endmjd = _qa.totime(m[1] + '/' + m[2])
                # date increment
                elif re.search(r'Step-size', line):
                    m = re.match(r'^[>\s]*Step-size\s+\S+\s+(\S+)\s+(\w+)', line)
                    if m:
                        unit = m[2]
                        if unit == 'minutes':
                            # this is the default unit returned by the JPL-Horizons
                            theunit = 'min'
                        elif unit == 'hours':
                            theunit = 'h'
                        elif unit == 'days':
                            theunit = 'd'
                        elif unit == 'steps':
                            theunit = 'steps'
                        elif unit == 'calendar':
                            raise RuntimeError('Unit of Step-size in calendar month or year is not supported.')
                        else:
                            raise RuntimeError(f'Unit of Step-size, {unit} is not recognized.')
                        if theunit == 'd':
                            dmjd = m[1]
                        elif theunit == 'steps':
                        #    print('endtime=',endmjd['value'])
                            dmjd = (endmjd['value'] - startmjd['value'])/int(m[1])
                        else:
                            dmjd = _qa.convert(_qa.totime(m[1] + theunit), 'd')['value']
                        headerdict['dMJD'] = dmjd
                        if startmjd is not None:  # start mjd should be available before step-size line
                            # MJD0 = firstMJD - dMJD (as defined casacore MeasComet documentation)
                            headerdict['MJD0'] = startmjd['value'] - dmjd
                elif re.search(r'Center geodetic', line):
                    m = re.match(r'^[>\s]*Center geodetic\s*: ([-+0-9.]+,\s*[-+0-9.]+,\s*[-+0-9.]+)', line)
                    if m:
                        long_lat_alt = m[1].split(',')
                        headerdict['GeoLong'] = float(long_lat_alt[0])
                        headerdict['GeoLat'] = float(long_lat_alt[1])
                        headerdict['GeoDist'] = float(long_lat_alt[2])
                        # obs location
                elif re.search(r'Center-site name', line):
                    m = re.match(r'^[>\s]*Center-site name:(\s*)([a-zA-Z\s*\/\-\(\)\,]+)', line)
                    if m:
                        # For the topocentric location, currently only ALMA, VLA and GBT are translated 
                        # to the proper observatory name which recongnized by Measures.
                        # The key part is the name used by JPL-Horizons.
                        observatories = {'ALMA':'ALMA', 'VLA':'VLA', 'Green Bank':'GBT'} 
                        if m[2]:
                           headerdict['obsloc'] = m[2] 
                           for obs in observatories:
                               if obs in m[2]:
                                   headerdict['obsloc'] = observatories[obs]
                                   break
                        # Assume the query is made in ICRF reference frame
                        # and for casacore measures, this will be 'ICRS'
                        headerdict['posrefsys'] =  'ICRS'
                elif re.search(r'Target radii', line):
                    m = re.match(r'^[>\s]*Target radii\s*:\s*([0-9.]+\s*x\s*[0-9.]+\s*x\s*[0-9.]+)\s*km.*|'
                                 '^[>/s]*Target radii\s*:\s*([0-9.]+)\s*km', line)
                    if m:
                        matchlist = m.groups()
                        radiiarr = np.zeros(3)
                        if len(matchlist)==2:
                            if m[2] is None:
                                radiiarr = np.asarray(np.array(m[1].split(' x ')), dtype=np.float64)
                                headerdict['radii'] = {'unit': 'km', 'value': radiiarr}
                            elif m[1] is None:
                                radiiarr = np.array([m[2],m[2],m[2]], dtype=np.float64)
                            headerdict['radii'] = {'unit': 'km', 'value': radiiarr}
                        else:
                            casalog.post(f"Unexpected number or matches for Target radii:{m.groups} (expected 2)", 'WARN')
                #rotational period (few pattens seem to exist)
                elif re.search(r'rot. period|Rotational period', line):
                    m = re.search(r'rot. period\s+\S*=\s*([0-9.]+h\s*[0-9.]+m\s*[0-9.]+\s*s)|'
                                    'rot. period\s+\S*=\s*([0-9.]+)(?:\s*\+-[0-9.]+)?\s*([dh])|'
                    'Rotational period\s*=\s+Synchronous', line)
                    if m:
                        if m[0].find('Synchronous') != -1:
                            headerdict['rot_per'] = 'Synchronous'
                        else:
                            if len(m.groups()) == 3:
                                if m[1] is None:
                                    headerdict['rot_per'] = _qa.quantity(m[2] + m[3])
                                elif m[2] is None and m[3] is None:
                                    #subm = re.search(r'([0-9]+)h\s*([0-9]+)m\s*([0-9.]+)\s*s',m[1])
                                    headerdict['rot_per'] = _qa.convert(re.sub(r'\s+','',m[1]),'h')
                # another variation of rotational period entry
                elif re.search(r'rot.\s+per.', line):
                    m = re.search(r'rot.\s+per.\s*=\s*([0-9.]+)\s*([dh])',line)
                    if m:
                        headerdict['rot_per'] = _qa.quantity(m[1]+m[2])
                # rotational period for asteroids
                elif re.search(r'ROTPER', line):
                    m = re.search(r'ROTPER=\s*([0-9.]+)\s*', line)
                    if m:
                        headerdict['rot_per'] = {'unit': 'h', 'value': float(m[1])}
                elif re.search(r'orbit period|orb per|orb. per.|orbital period', line.lower()) \
                        and not foundorbper:
                    m = re.search(r'Orbital period\s*[=~]\s*([-0-9.]+)\s*(\w+)\s*|'
                                  'orb. per., (\w)\s+=\s+([0-9.]+)\s+|'
                                  'orb per\s+=\s+([0-9.]+)\s+(\w+)\s+|'
                                  'orbit period\s+=\s+([0-9.]+)\s+(\w+)\s+', line)
                    if m:
                        if m[0].find('Orbital period') != -1:
                            headerdict['orb_per'] = {'unit': m[2], 'value': float(m[1])}
                        elif m[0].find('orb. per.') != -1:
                            headerdict['orb_per'] = {'unit': m[3], 'value': float(m[4])}
                        elif m[0].find('orb per') != -1:
                            headerdict['orb_per'] = {'unit': m[6], 'value': float(m[5])}
                        elif m[0].find('orbit period') != -1:
                            headerdict['orb_per'] = {'unit': m[8], 'value': float(m[7])}
                        if 'orb_per' in headerdict:
                            foundorbper = True
                elif re.search(r'Mean [Tt]emperature', line):
                    m = re.search(r'Mean [Tt]emperature\s+\(K\)\s+=\s+([0-9.]+)\s+',line)
                    if m:
                        headerdict['T_mean'] = {'unit':'K', 'value':float(m[1])}
                # start reading data
                elif re.search(r'\s*Date__\(UT\)', line):
                    incolnames = line.split()
                elif re.search(r'\$\$SOE', line):
                    readthedata = True
                elif re.search(r'\$\$EOE', line):
                    readthedata = False
                elif readthedata:
                    datalines += 1
                    outfile.write(line + '\n')
                elif re.search(r'Projected output length', line):
                    # this message occurs when requested data is too large
                    exceedthelinelimit = line
                elif re.search(r'Multiple major-bodies match', line):
                    ambiguousname = line
                else:
                    pass
                lcnt += 1
            if 'radii' in headerdict:
                radiival = headerdict['radii']['value']
                meanrad = _mean_radius(radiival[0], radiival[1], radiival[2])
                headerdict['meanrad'] = {'unit': 'km', 'value': meanrad}
            if datalines == 0:
                casalog.post("No ephemeris data was found", "WARN")
                queryerrmsg = exceedthelinelimit if exceedthelinelimit != None else ambiguousname
                
                if queryerrmsg:
                    raise RuntimeError(f'Error occurred at the query:{queryerrmsg}')
                else:
                    raise RuntimeError(f'Error occurred at the query')
            else:
                casalog.post(f'Number of data lines={datalines}')
            casalog.post(f'Number of all lines in the file={lcnt}')
            #print("headerdict=", headerdict)
        # output to a casa table

        # check the input data columns and stored the order as indices
        foundncols = 0
        indexoffset = 0
        colkeys = {}
        #print('incolnames=',incolnames)
        if incolnames is not None:
            for outcolname in cols:
                # all colnames in cols should have unit defined.
                if 'unit' in cols[outcolname]:
                    colkeys[outcolname] = np.array([cols[outcolname]['unit']])
                else:
                    raise KeyError(f'Missing unit for {outcolname}.')
                inheadername = cols[outcolname]['header']
                # expect date is in the first column (date and mm:hh seperated by spaces)
                if outcolname == 'MJD':
                    for incol in incolnames:
                        if re.search(inheadername + '.+', incol):
                            cols[outcolname]['index'] = incolnames.index(incol)
                            foundncols += 1
                            indexoffset = 1
                    if 'index' not in cols[outcolname]:
                        casalog.post("Cannot find the Date column", "WARN")
                elif outcolname == 'RA':
                    for incol in incolnames:
                        if re.search(inheadername + '.+(ICRF).+', incol):
                            cols[outcolname]['index'] = incolnames.index(incol) + indexoffset
                            foundncols += 1

                    if 'index' not in cols[outcolname]:
                        casalog.post("Cannot find the astrometric RA and Dec column", "WARN")
                elif outcolname == 'DEC':
                    if 'index' in cols['RA']:
                        # Dec data col is next to RA data col
                        cols[outcolname]['index'] = cols['RA']['index'] + 1
                        # add additional offset for index (4 data columns at this point)
                        indexoffset +=1
                        foundncols += 1
                else:
                    if inheadername in incolnames:
                        cols[outcolname]['index'] = incolnames.index(inheadername) + indexoffset
                        foundncols += 1
                    else:
                        casalog.post(f"Cannot find {inheadername}", "WARN")

            #print(cols)
            casalog.post(f"expected n cols = {len(cols)} ")
            casalog.post(f"foundncols = {foundncols}")
            if foundncols == len(cols):
                # Format the data to comply with measure/setjy
                casalog.post("Found all the required columns")
                with open(tempconvfname, 'w') as outf, open(tempfname, 'r') as inf:
                    ndata = 0
                    earliestmjd = None
                    mjd = None
                    prevmjd = None
                    for line in inf:
                        outline = ''
                        sep = ' '
                        extraindexoffset = 0 
                        rawtempdata = line.split()
                        tempdata = ['-999.0' if x == 'n.a.' else x for x in rawtempdata]
                        # Check if there are extra codes (solar/lunar presence codes) after the date.
                        # The extra codes may be present for topocentric observer location. Even in that case
                        # not every data line has the code so it needs to be checked line by line.
                        ntempdata = len(tempdata)
                        nincolnames = len(incolnames) + 2 # 2 for date_time and ra_dec splits 
                        if ntempdata > nincolnames :
                            # assume there is the extra codes (no space between solar and lunar presence code)
                            if ntempdata - nincolnames == 1:
                                extraindexoffset = 1
                            else:
                                raise RuntimeError(f'Unexpected number of data items was detected:'+
                                                    f' {ntempdata} while expected {nincolnames}') 
                        # construct mjd from col 1 (calendar date) + col 2 (time)
                        caldatestr = tempdata[0] + ' ' + tempdata[1]
                        mjd = _qa.totime(caldatestr)
                        if mjd == prevmjd:
                           raise RuntimeError(f'Duplicated timestamp, {mjd}, is detected. This may occur when '+
                            'a time range is specified in MJD with a short time interval. If that is the case, '+
                            'try calendar date+time string for the time range.') 
                        outline += str(mjd['value']) + sep
                        # position
                        rad = tempdata[cols['RA']['index'] + extraindexoffset]
                        decd = tempdata[cols['DEC']['index'] + extraindexoffset]
                        outline += rad + sep + decd + sep
                        # geocentric dist. (Rho)
                        delta = tempdata[cols['Rho']['index'] + extraindexoffset]
                        outline += delta + sep
                        # geocentric range rate (RadVel)
                        valinkmps = tempdata[cols['RadVel']['index'] + extraindexoffset]
                        deldot = _qa.convert(_qa.quantity(valinkmps+'km/s'), 'AU/d' )['value']
                        outline += str(deldot) + sep
                        # NP_ang & NP_dist
                        npang = tempdata[cols['NP_ang']['index'] + extraindexoffset]
                        npdist = tempdata[cols['NP_dist']['index'] + extraindexoffset]
                        outline += npang + sep + npdist + sep
                        # DiskLong & DiskLat
                        disklong = tempdata[cols['DiskLong']['index'] + extraindexoffset]
                        disklat = tempdata[cols['DiskLat']['index'] + extraindexoffset]
                        outline += disklong + sep + disklat + sep
                        # sub-long & sub-lat
                        sllon = tempdata[cols['Sl_lon']['index'] + extraindexoffset]
                        sllat = tempdata[cols['Sl_lat']['index'] + extraindexoffset]
                        outline += sllon + sep + sllat + sep
                        # r, rot
                        r = tempdata[cols['r']['index'] + extraindexoffset]
                        rdot = tempdata[cols['rdot']['index'] + extraindexoffset]
                        outline += r + sep + rdot + sep
                        # S-T-O
                        phang = tempdata[cols['phang']['index'] + extraindexoffset]
                        outline += phang
                        outf.write(outline + '\n')

                        ndata += 1
                        if ndata == 1:
                            earliestmjd = mjd
                        prevmjd = mjd
                    # record first and last mjd in the data
                    headerdict['earliest'] = _me.epoch('UTC',earliestmjd)
                    headerdict['latest'] = _me.epoch('UTC', mjd)

            else:
                missingcols = len(cols) - foundncols
                casalog.post(r"Missing {missingcols}")

            # final step: convert to a CASA table
            dtypes = np.array(['D' for _ in range(len(cols))])
            _tb.fromascii(outtable, tempconvfname, sep=' ', columnnames=list(cols.keys()),
                          datatypes=dtypes.tolist())
            _tb.done()
            
            # fill keyword values in the ephem table
            if os.path.exists(outtable):
                _fill_keywords_from_dict(headerdict, colkeys, outtable)
                casalog.post(f"Output is written to a CASA table, {outtable}")
            else:
                raise RuntimeError("Error occured. The output table, " + outtable + "is not generated")
        else: #incolnames is None
            raise RuntimeError("No data header was found.")

    except RuntimeError as e:
        casalog.post(str(e),"SEVERE")
    
    finally:
        tempfiles = [tempfname, tempconvfname]
        _clean_up(tempfiles)

# This needs to change to read the ephem data from the generated casa table.
# Will be called in fill_keywords_from_dict()

# Copied from JPLephem_reader2.py for mean_radius calc (when no NP-ang, NP_dist)
# Note: this is not fully verified for the correctness but has been used 
# through JPLephem_reader2 for processing the ephemeris data to a CASA table.
# For setjy, solar_system_setjy does own mean radius calculation and 
# meanrad is the table is not used for that. 2022.01.12 TT
def _mean_radius(a, b, c):
    """
    Return the average apparent mean radius of an ellipsoid with semiaxes
    a >= b >= c.
    "average" means average over time naively assuming the pole orientation
    is uniformly distributed over the whole sphere, and "apparent mean radius"
    means a radius that would give the same area as the apparent disk.
    """
    # This is an approximation, but it's not bad.
    # The exact equations for going from a, b, c, and the Euler angles to the
    # apparent ellipse are given in Drummond et al., Icarus, 1985a.
    # It's the integral over the spin phase that I have approximated, so the
    # approximation is exact for b == a and appears to hold well for b << a.
    R = 0.5 * c ** 2 * (1.0 / b ** 2 + 1.0 / a ** 2)  # The magic ratio.
    if R < 0.95:
        sqrt1mR = sqrt(1.0 - R)
        # There is fake singularity (RlnR) at R = 0, but it is unlikely to
        # be a problem.
        try:
            Rterm = 0.5 * R * log((1.0 + sqrt1mR) / (1.0 - sqrt1mR)) / sqrt1mR
        except RuntimeError:
            Rterm = 0.0
    else:
        # Use a (rapidly converging) series expansion to avoid a fake
        # singularity at R = 1.
        Rterm = 1.0  # 0th order
        onemR = 1.0 - R
        onemRtothei = 1.0
        for i in range(1, 5):  # Start series at 1st order.
            onemRtothei *= onemR
            Rterm -= onemRtothei / (0.5 + 2.0 * i ** 2)
    avalfabeta = 0.5 * a * b * (1.0 + Rterm)
    return sqrt(avalfabeta)


def _mean_radius_with_known_theta(disklat, radii):
    """
    Mean radius calculation extracted from solar_system_setjy
    """
    Req = 1000.0 * (radii[0] + radii[1]) / 2
    Rp = 1000.0 * radii[2]

    Rmean = 0.0
    index = 0
    for index, lat in enumerate(disklat):
        if lat == -999.0:
            lat = 0.0
        rlat = lat * (180. / np.pi)
        Rpap = sqrt(Req * Req * sin(rlat) ** 2.0 +
                    Rp * Rp * cos(rlat) ** 2.0)
        #    Rmean += ( sqrt (Rpap * Req) - Rmean )/ (index + 1.0)
        Rmean += sqrt(Rpap * Req)

    return (Rmean / (index + 1.0)) / 1000.0


def _fill_keywords_from_dict(keydict, colkeys, tablename):
    # open tb
    # get ra,dec,np_ra, np_dec, and radii in the keyword
    # call mod version of mean_radius_with_known_theta
    orderedmainkeys = ['VS_CREATE','VS_DATE','VS_TYPE','VS_VERSION','NAME',
                       'MJD0','dMJD','GeoDist','GeoLat','GeoLong','obsloc',
                       'posrefsys','ephemeris_source','earliest','latest','radii',
                       'meanrad','orb_per','rot_per','T_mean']
    try:
        _tb.open(tablename, nomodify=False)
        disklat = _tb.getcol('DiskLat')
        if 'radii' in keydict:
            radiival = keydict['radii']['value']
            calc_meanrad = _mean_radius_with_known_theta(disklat, radiival)
            if calc_meanrad and ('meanrad' in keydict):
                keydict['meanrad']['value'] = calc_meanrad
        if 'rot_per' in keydict and 'orb_per' in keydict and keydict['rot_per'] == 'Synchronous':
            keydict['rot_per'] = keydict['orb_per']
        sortedkeydict = {k: keydict[k] for k in orderedmainkeys if k in keydict}
        #print('sorteddict=',sortedkeydict)
        for k in sortedkeydict:
            _tb.putkeyword(k, sortedkeydict[k])
        datacolnames = _tb.colnames()
        for col in datacolnames:
            if col in colkeys:
                _tb.putcolkeyword(col, 'QuantumUnits', colkeys[col])
        # add table info required by measComet
        maintbinfo = {'readme': 'Derivied by jplhorizons-query from JPL-Horizons API '
                                '(https://ssd.jpl.nasa.gov/api/horizons.api)',
                      'subType':'Comet',
                      'type': 'IERS'}
        _tb.putinfo(maintbinfo)
        _tb.flush()
        _tb.done()
    except RuntimeError:
        casalog.post('Internal error: Cannot add the data in keywords', "WARN")

def _clean_up(filelist):
    """
    Clean up the temporary files 
    """
    for f in filelist:
        if os.path.exists(f):
            os.remove(f)
