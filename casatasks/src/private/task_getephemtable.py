import os
import re
import shutil
from casatasks import casalog
from casatasks.private.jplhorizons_query import gethorizonsephem


def getephemtable(objectname, asis, timerange, interval, outfile, rawdatafile, overwrite):
    """Retrieve the ephemeris data of a specific ephemeris object by sending
    a query to JPL's Horizons system and creates the ephemeris data stored in a CASA table format.
    """
    casalog.origin('getephemtable')

#Python script

    if isinstance(objectname, str):
        if not objectname.strip():
            raise ValueError("objectname must be specified")

        match = re.match(r'(\s*)([0-9]+)(\s*)', objectname)
        if match is not None and not asis:
            raise RuntimeError("objectname is given as an ID number, need to set asis=True")

    # split timerange to start and end times
    if isinstance(timerange, str):
        if not timerange.strip():
            raise ValueError("timerange must be specified")

        if timerange.find('~') != -1:
            timerange = timerange.replace(' ','')
            (starttime, stoptime) = timerange.split('~')
            starttime = starttime.upper()
            stoptime = stoptime.upper()
            if starttime.startswith('JD'):
                if not stoptime.startswith('JD'):
                    try:
                        float(stoptime)
                        stoptime = 'JD' + stoptime
                    except Exception as e:
                        raise ValueError("Error translating stop time of timerange specified in JD.") from e
            # JPL-Horizons does not accept MJD for time specification.
            elif starttime.startswith('MJD'):
                try:
                    starttime = 'JD'+str(float(starttime.strip('MJD')) + 2400000.5)
                except Exception as e:
                    raise ValueError("Error translating start time of timerange specified in MJD.") from e
                if not stoptime.startswith('JD'):
                    if stoptime.startswith('MJD'):
                        stoptime = stoptime.strip('MJD')
                    try:
                        stoptime = 'JD' + str(float(stoptime) + 2400000.5)
                    except Exception as e:
                        raise ValueError("Error translating stop time of timerange specified in MJD.") from e
            else:
                matchstart = re.match(r'(\s*)([0-9][0-9][0-9][0-9])\/([0-9][0-9])\/([0-9][0-9])([/:0-9]*)',starttime)
                matchstop = re.match(r'(\s*)([0-9][0-9][0-9][0-9])\/([0-9][0-9])\/([0-9][0-9])([/:0-9]*)',stoptime)
                #print(f'startime={starttime}, stoptime={stoptime}, matchstart={matchstart}, matchstop={matchstop}')
                if matchstart is None or matchstop is None:
                    raise ValueError("Error in timerange format. Use YYYY/MM/DD/hh:mm or Julian date with a prefix 'JD' Modified Julian date with a prefix 'MJD'")
        else:
            raise ValueError("timerange needs to be specified with starttime and stoptime delimited by ~ .")


    # check for interval
    if isinstance(interval, str):
        match = re.match(r'(\d+(?:\.\d*)?)(\s*)([a-zA-Z]*)', interval)
        if match is not None:
            (valstr, _, unitstr) = list(match.groups())
            if valstr.isdigit():
                if unitstr != '':
                    intervalstr = valstr+unitstr
                else:
                    intervalstr = valstr
            else:
                raise ValueError("interval value must be integer")
        else:
            raise ValueError("interval must be a string with integer with unit ")

    # outfile and rawdatafile check
    if not outfile.strip():
        raise ValueError("outfile must be specified")
    elif os.path.exists(outfile):
        if not overwrite:
            raise ValueError(f'{outfile} exists and overwrite=False')
        else:
            casalog.post(f'{outfile} exists, will be overwritten', 'WARN')
        shutil.rmtree(outfile)
    if os.path.exists(rawdatafile):
        if not overwrite:
            raise ValueError(f'{rawdatafile} exists and overwrite=False')
        else:
            casalog.post(f'{rawdatafile} exists, will be overwritten', 'WARN')

    # call the JPL-Horizons query function
    gethorizonsephem(objectname, starttime, stoptime, intervalstr, outfile, asis, rawdatafile)
