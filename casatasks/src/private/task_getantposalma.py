from casatasks import casalog
from casatools import quanta
import certifi
from datetime import datetime
import json, os, shutil
import ssl
from urllib import request
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, urlparse


def _is_valid_url_host(url):
    parsed = urlparse(url)
    return bool(parsed.netloc)


def _query(url):
    myjson = None
    response = None
    try:
        context = ssl.create_default_context(cafile=certifi.where())
        with request.urlopen(url, context=context, timeout=400) as response:
            if response.status == 200:
                myjson = response.read().decode('utf-8')
    except HTTPError as e:
        casalog.post(
            f"Caught HTTPError: {e.code} {e.reason}: {e.read().decode('utf-8')}",
            "WARN"
        )
    except URLError as e:
        casalog.post(f"Caught URLError: {str(e)}", "WARN")
    except Exception as e:
        casalog.post(f"Caught Exception when trying to connect: {str(e)}", "WARN")
    return myjson


def getantposalma(
    outfile='', overwrite=False, asdm='', tw='', snr="default", search='both_latest',
    hosts=['tbd1.alma.cl', 'tbd2.alma.cl']
):
    r"""
Retrieve antenna positions by querying ALMA web service.

[`Description`_] [`Examples`_] [`Development`_] [`Details`_]


Parameters
   - outfile_ (path='') - Name of output file to which to write retrieved antenna positions.
   - overwrite_ (bool=False) - Overwrite a file by the same name if it exists?
   - asdm_ (string='') - The associated ASDM name. Must be specified
   - tw_ (string='') - Optional time window in which to consider baseline measurements in the database, when calculating the antenna positions.
   - snr_ (float=0) - Optional signal-to-noise.
   - search_ (string='both_latest') - Search algorithm to use.
   - hosts_ (stringVec=['https://asa.alma.cl/uncertainties-service/uncertainties/versions/last/measurements/casa/']) - Priority-ranked list of hosts to query.




.. _Description:

Description

.. warning:: **WARNING**: This task should be considered experimental
   since the values returned by the JAO service are in the process of
   being validated.

This task retrieves ALMA antenna positions via a web service which runs
on an ALMA-hosted server. The antenna positions are with respect to ITRF.
The user must specify the value of the outfile parameter. This parameter
is the name of the file to which the antenna positions will be written.
This file can then be read by gencal so that it can use the most up to
date antenna positions for the observation.

The web service is described by the server development team and can be
found `at this location <https://asw.alma.cl/groups/ASW/-/packages/843>`__. 

The input parameters are discussed in detail below.

outfile is required to be specified. It is the name of the file to which to
write antenna positions.

overwrite If False and a file with the same name exists, and exception
will be thrown. If true, an existing file with the same name will be
overwriten.

asdm is required to be specified. It is the associated ASDM UID in the
form uid://A002/Xc02418/X29c8. 

tw is an optional parameter. It is time window in which the antenna positions
are required, specified as a comma separated pair. Times are UTC and are
expressed in YY-MM-DDThh:mm:ss.sss format. The end time must be later than
the begin time.

snr is an optional parameter. If changed from the default value "default", 
it must be a nonnegative number representing the signal-to-noise ratio. Antenna
positions which have corrections less than this value will not be written.
If not specified, the default snr as defined by the web service will be used.
The server side default value may change over time as determined by the server
side (non-CASA) team. As of this writing (March 2025), the web service team has
not provided publicly facing documentation on the details of how the default
value is chosen. The most recent information they have provided to us is that
the default value is 5.0.

tw and search are optional parameters and are coupled as follows. search
indicates the search algorithm to use to find the desired antenna positions.
Supported values of this parameter at the time of writing are 'both_latest'
and 'both_closest'. The task passes the value of the search parameter verbatim to
the web service, meaning that users can take advantage of new search algorithms
as the web service team brings them online. The default algorithm used is
'both_latest'. In general, the search is limited in time to the specified
value of tw (time window). However, in the case that tw is not specified, the
following rules apply. For 'both_latest', the last updated position for each
antenna within the specified time window, or, if tw is not specified, within
30 days after the observation will be returned, taking into account snr if
specified, if provided.

For 'both_closest', if tw is not specified, the position
of each antenna closest in time to the observation, within 30 days (before
or after the observation) will be returned, subject to the value of snr if it
is specified. 

hosts is a required parameter. It is a list of hosts to query, in order of
priority, to obtain positions. The first server to respond with a valid result is
the only one that is used. That response will be written and no additional
hosts will be queried.

The format of the returned file is a two element list encoded in json. The first
element is a stringfied dictionary that contains antenna names as keys, with each
value being a three element list of x, y, and z coordinates in ITRF. The second
element is a dictionary containing various (possibly helpful) metadata that were
used when the task was run. The following code may be used to load these data
structures into python variables.
    
    ::
        
        import ast, json
        ...
        with open("outfile.json", "r") as f:
            antpos_str, md_dict = json.load(f)
            antpos_dict = ast.literal_eval(antpos_str)


.. _Examples:

Examples
   Get antenna positions which have positions with a signal-to-noise ratio
   greater than 5.
   
   ::
   
      getantposalma(
          outfile='my_ant_pos.json', asdm='valid ASDM name here', snr=5,
          hosts=['tbd1.alma.cl', 'tbd2.alma.cl']
     )
   

.. _Development:

Development
   No additional development details




.. _Details:


Parameter Details
   Detailed descriptions of each function parameter

.. _outfile:

| ``outfile (path='')`` - Name of output file to which to write antenna positions. If a file by this name already exists, it will be silently overwritten. The written file will be in JSON format.
|    default: none
|    Example: outfile='my_alma_antenna_positions.json'

.. _overwrite:

| ``overwrite (bool=False)`` - Overwrite a file by the same name if it exists? If False and a file
|                with the same name exists, and exception will be thrown.

.. _asdm:

| ``asdm (string='')`` - The associated ASDM name. Must be specified. The ASDM is not required to be on the file system; its value is simply passed to the web service.
|                       default: ''
|                       Example:asdm='uid://A002/X10ac6bc/X896d'

.. _tw:

| ``tw (string='')`` - Optional time window in which to consider baseline measurements in the database, when calculating the antenna positions. Format is of the form begin_time,end_time, where times must be specified in YYYY-MM-DDThh:mm:ss.sss format and end_time must be later than begin time. Times should be UTC.
|          Example: tw='2023-03-14T00:40:20,2023-03-20T17:58:20'

.. _snr:

| ``snr (float=0)`` - Optional signal-to-noise. Antenna positions which have corrections with S/N less than this value will not be retrieved nor written. If not specified, positions of all antennas will be written.
|            default: 0 (no snr constraint will be used) 
|            Example: snr=5.0

.. _search:

| ``search (string='both_latest')`` - Search algorithm to use. Supported values are "both_latest" and "both_closest". For "both_latest", the last updated position for each antenna within 30 days after the observation will be returned, taking into account snr if specified. If provided, tw will override the 30 day default value. For "both_closest", the position of each antenna closest in time to the observation, within 30 days (before or after the observation) will be returned, subject to the value of snr if it is specified. If specified, the value of tw will override the default 30 days. The default algorithm to use will be "both_latest".
|          Example: search="both_closest"

.. _hosts:

| ``hosts (stringVec=['https://asa.alma.cl/uncertainties-service/uncertainties/versions/last/measurements/casa/'])`` - Priority-ranked list of hosts to query to obtain positions. Only one server that returns a list of antenna positions is required. That response will be written and no additional hosts will be queried.
|            Example: hosts=["server1.alma.cl", "server2.alma.cl"]


    """
    if not outfile:
        raise ValueError("Parameter outfile must be specified")
    md = {
        "caltype": "ALMA antenna positions",
        "description": "ALMA ITRF antenna positions in meters",
        "product_code": "antposalma",
        "outfile": outfile
    }
    if not overwrite and os.path.exists(outfile):
        raise RuntimeError(
            f"A file or directory named {outfile} already exists and overwrite "
            "is False, so exiting. Either rename the existing file or directory, "
            "change the value of overwrite to True, or both."
        )
    if not hosts:
        raise ValueError("Parameter hosts must be specified")
    if isinstance(hosts, list) and not hosts[0]:
        raise ValueError("The first element of the hosts list must be specified")
    md["hosts"] = hosts
    _qa = quanta()
    parms = {}
    if asdm:
        parms['asdm'] = asdm
    else:
        raise ValueError("parameter asdm must be specified")
    if tw:
        z = tw.split(",")
        if len(z) != 2:
            raise ValueError(
                "Parameter tw should contain exactly one comma that separates two times"
            )
        s0, s1 = z
        msg = "The correct format is of the form YYYY-MM-DDThh:mm:ss."
        try:
            t_start = _qa.quantity(_qa.time(s0, form="fits")[0])
        except Exception as e:
            raise ValueError(f"Begin time {s0} does not appear to have a valid format. {msg}")
        try:
            t_end = _qa.quantity(_qa.time(s1, form="fits")[0])
        except Exception as e:
            raise ValueError(f"End time {s1} does not appear to have a valid format. {msg}")
        if _qa.ge(t_start, t_end):
            raise ValueError(
                f"Parameter tw, start time ({z[0]}) must be less than end time ({z[1]})."
            )
        parms["tw"] = tw
    if isinstance(snr, str):
        if snr != "default":
            raise ValueError("If snr is a string, it's only permissible value is 'default'")
    elif snr < 0.0:
        raise ValueError(f"If a number, parameter snr ({snr}) must be non-negative.")
    elif snr >= 0.0:
        parms["snr"] = snr
    if search:
        parms['search'] = search
    qs = f"?{urlencode(parms)}"
    md.update(parms)
    antpos = None
    for h in hosts:
        if not _is_valid_url_host(h):
            raise ValueError(
                f'Parameter hosts: {h} is not a valid host expressed as a URL.'
            )
        url = f"{h}/{qs}"
        casalog.post(f"Trying {url} ...", "NORMAL")
        antpos = _query(url)
        if antpos:
            md["successful_url"] = url
            antpos = json.loads(antpos)
            break
    if not antpos:
        raise RuntimeError("All URLs failed to return an antenna position list.")
    if os.path.exists(outfile):
        if overwrite:
            if os.path.isdir(outfile):
                casalog.post(
                    f"Removing existing directory {outfile} before writing new "
                    "file of same name",
                    "WARN"
                )
                shutil.rmtree(outfile)
            else:
                casalog.post(
                    f"Removing existing file {outfile} before writing now file of "
                    "same name",
                    "WARN"
                )
                os.remove(outfile)
        else:
            raise RuntimeError(
                "Logic Error: shouldn't have gotten to this point with overwrite=False"
            )
    md["timestamp"] = str(datetime.now())
    with open(outfile, "w") as f:
        json.dump({"data": antpos, "metadata": md}, f)
