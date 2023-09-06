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
       # split timerange to start and end times
       if type(timerange) == str:
           if timerange.find('~'):
               timerange = timerange.replace(' ','')
               print('timerange before split =', timerange)
               (starttime, stoptime) = timerange.split('~')
               starttime = starttime.upper()
               stoptime = stoptime.upper() 
               if starttime.startswith('JD'):
                   if not stoptime.startswith('JD'):
                       try:
                           float(stoptime)
                           stoptime = 'JD' + stoptime
                       except:
                           raise TypeError("Error translating stop time of timerange specified in JD.")
               # JPL-Horizons does not accept MJD for time specification.
               elif starttime.startswith('MJD'):
                  try:
                      starttime = 'JD'+str(float(starttime.strip('MJD')) + 2400000.5) 
                  except:
                      raise TypeError("Error translating start time of timerange specified in MJD.")
                  print("stoptime=",stoptime)
                  if not stoptime.startswith('JD'):
                      if stoptime.startswith('MJD'):
                          stoptime = stoptime.strip('MJD')
                      try:
                          stoptime = 'JD' + str(float(stoptime) + 2400000.5)
                      except:
                          raise TypeError("Error translating stop time of timerange specified in MJD.")
           else:
               raise TypeError("timerange needs to be specified with starttime and stoptime connected by ~ .")
 

       # check for interval
       if type(interval) == str:
          match = re.match(r'([0-9]+)(\s*)([a-zA-Z]+)', interval)
          if match is not None:
             (intstr, ws, unitstr) = list(match.groups())
             intervalstr = intstr+unitstr
          else:
             raise TypeError("interval must contains integer value and unit")

       # outfile and rawdatafile check
       if os.path.exists(outfile):
            casalog.post(f'{outfile} exists, will be overwritten', 'WARN')
       if os.path.exists(rawdatafile):
            casalog.post(f'{rawdatafile} exists, will be overwritten', 'WARN')

       # call the JPL-Horizons query function
       gethorizonsephem(objectname, starttime, stoptime, intervalstr, outfile, asis, rawdatafile)

