from __future__ import absolute_import
import os
import re
from casatasks import casalog
from .jplhorizons_query import gethorizonsephem


def getephemtable(objectname, asis, timerange, interval, outfile, rawdatafile):
       """Retrieve the ephemeris data of a specific ephemeris object by sending
       a query to JPL's Horizons system and creates the ephemeris data stored in a CASA table format.
       """
       casalog.origin('getephemtable')

       #Python script

       if type(objectname) == str: 
          if not objectname.strip():
              raise ValueError("objectname must be specified")

          match = re.match(r'(\s*)([0-9]+)(\s*)', objectname)
          if match is not None and not asis:
              raise RuntimeError("objectname is given as an ID number, need to set asis=True")

       # split timerange to start and end times
       if type(timerange) == str:
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
                       except:
                           raise ValueError("Error translating stop time of timerange specified in JD.")
               # JPL-Horizons does not accept MJD for time specification.
               elif starttime.startswith('MJD'):
                  try:
                      starttime = 'JD'+str(float(starttime.strip('MJD')) + 2400000.5) 
                  except:
                      raise ValueError("Error translating start time of timerange specified in MJD.")
                  print("stoptime=",stoptime)
                  if not stoptime.startswith('JD'):
                      if stoptime.startswith('MJD'):
                          stoptime = stoptime.strip('MJD')
                      try:
                          stoptime = 'JD' + str(float(stoptime) + 2400000.5)
                      except:
                          raise ValueError("Error translating stop time of timerange specified in MJD.")
               else:
                  matchstart = re.match(r'(\s*)([0-9][0-9][0-9][0-9])\/([0-9][0-9])\/([0-9][0-9])([/:0-9]*)',starttime)
                  matchstop = re.match(r'(\s*)([0-9][0-9][0-9][0-9])\/([0-9][0-9])\/([0-9][0-9])([/:0-9]*)',stoptime)
                  print(f'startime={starttime}, stoptime={stoptime}, matchstart={matchstart}, matchstop={matchstop}')
                  if matchstart is None or matchstop is None:
                      raise ValueError("Error in timerange format. Use YYYY/MM/DD/hh:mm or Julian date with a prefix 'JD' Modified Julian date with a prefix 'MJD'")
           else:
               raise ValueError("timerange needs to be specified with starttime and stoptime connected by ~ .")
 

       # check for interval
       if type(interval) == str:
          match = re.match(r'([0-9]+)(\s*)([a-zA-Z]+)', interval)
          if match is not None:
             (intstr, ws, unitstr) = list(match.groups())
             intervalstr = intstr+unitstr
          else:
             raise ValueError("interval must contains integer value and unit")

       # outfile and rawdatafile check
       if not outfile.strip():
           raise ValueError("outfile must be specified")
       elif os.path.exists(outfile):
            casalog.post(f'{outfile} exists, will be overwritten', 'WARN')

       if os.path.exists(rawdatafile):
            casalog.post(f'{rawdatafile} exists, will be overwritten', 'WARN')

       # call the JPL-Horizons query function
       gethorizonsephem(objectname, starttime, stoptime, intervalstr, outfile, asis, rawdatafile)

