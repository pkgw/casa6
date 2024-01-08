from casatasks import casalog
from casatools import quanta
import json, os, shutil
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
        with request.urlopen(url) as response:
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
    outfile='', overwrite=False, asdm='', tw='', snr=0, search='both_latest',
    hosts=['tbd1.alma.cl', 'tbd2.alma.cl']
):
    r"""
Retrieve antenna positions by querying ALMA web service.

[`Description`_] [`Examples`_] [`Development`_] [`Details`_]


Parameters
   - outfile_ (path='') - Name of output file to which to write retrieved antenna positions.
   - asdm_ (string='') - The associated ASDM name. Must be specified
   - tw_ (string='') - Optional time window to which to limit search for antenna positions.
   - snr_ (float=0) - Optional signal-to-noise.
   - search_ (string='both_latest') - Search algorithm to use.
   - hosts_ (stringVec=['tbd1.alma.cl', 'tbd2.alma.cl']) - Priority-ranked list of hosts to query.




.. _Description:

Description
Antpos retrieves ALMA antenna positions via a web service which runs on an
ALMA-hosted server. The antenna positions are with respect to ITRF. The
user must specify the value of the outfile parameter. This parameter is
the name of the file to which the antenna positions will be written. This
file can then be read by gencal so that it can use the most up to date
antenna positions for the observation.

The input parameters are discussed in detail below.

outfile is required to be specified. It is the name of the file to which to
write antenna positions. If a file with the same name exists it will be
silently overwritten.

asdm is required to be specified. It is the associated ASDM name.
tw is an optional parameter. It is time window in which the antenna positions
are required, specified as a comma separated pair. Times are UTC and are
expressed in YY-MM-DDThh:mm:ss.sss format. The end time must be later than
the begin time.

snr is an optional parameter. It is the signal-to-noise ratio. Antenna
positions which have corrections less than this value will not be written.
If not specified, positions of all antennas will be written.

search is an optional parameter, It is the search algorithm to use.
Supported values are 'both_latest' and 'both_closest'. For 'both_latest',
the last updated position for each antenna within 30 days after the
observation will be returned, taking into account snr if specified. If provided,

tw will override the 30 day default value. For 'both_closest', the position
of each antenna closest in time to the observation, within 30 days (before
or after the observation will be returned, subject to the value of snr if it
is specified. If specified, the value of tw will override the default 30 days.
The default algorithm used is 'both_latest'.

hostss is a required parameter. It is a list of hosts to query, in order of
priority, to obtain positions. The first server to respond with a valid result is
the only one that is used. That response will be written and no additional
hostss will be queried.


.. _Examples:

Examples
   Get antenna positions which have positions with a signal-to-noise ratio
   greater than 5.
   
   ::
   
      antpos(
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

.. _asdm:

| ``asdm (string='')`` - The associated ASDM name. Must be specified
|                       default: ''
|                       Example:asdm='uid___A002_X10ac6bc_X896d

.. _tw:

| ``tw (string='')`` - Optional time window to which to limit search for antenna positions. Format is of the form begin_time,end_time, where times must be specified in YYYY-MM-DDThh:mm:ss.sss format and end_time must be later than begin time. Times should be UTC.
|          Example: tw='2023-03-14T00:40:20,2023-03-20T17:58:20'

.. _snr:

| ``snr (float=0)`` - Optional signal-to-noise. Antenna positions which have corrections with S/N less than this value will not be retrieved nor written. If not specified, positions of all antennas will be written.
|            default: 0 (no snr constraint will be used) 
|            Example: snr=5.0

.. _search:

| ``search (string='both_latest')`` - Search algorithm to use. Supported values are "both_latest" and "both_closest". For "both_latest", the last updated position for each antenna within 30 days after the observation will be returned, taking into account snr if specified. If provided, tw will override the 30 day default value. For "both_closest", the position of each antenna closest in time to the observation, within 30 days (before or after the observation) will be returned, subject to the value of snr if it is specified. If specified, the value of tw will override the default 30 days. The default algorithm to use will be "both_latest".
|          Example: search="both_closest"

.. _hosts:

| ``hosts (stringVec=['tbd1.alma.cl', 'tbd2.alma.cl'])`` - Priority-ranked list of hosts to query to obtain positions. Only one server that returns a list of antenna positions is required. That response will be written and no additional hosts will be queried.
|            Example: hosts=["server1.alma.cl", "server2.alma.cl"]


    """
    if not outfile:
        raise ValueError("Parameter outfile must be specified")
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
    if snr < 0:
        raise ValueError(f"Parameter snr ({snr}) must be non-negative.")
    elif snr > 0:
        parms["snr"] = snr
    if search:
        if search in ["both_latest", "both_closest"]:
            parms["search"] = search
        else:
            raise ValueError(
                f"Parameter search (={search}) must have a value of either "
                "'both_latest' or 'both_closest'."
            )
    qs = f"?{urlencode(parms)}"
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
            break
    if not antpos:
        raise RuntimeError("All URLs failed to return an antenna position list.")
    if os.path.exists(outfile):
        if overwrite:
            if os.path.isdir(outfile):
                casalog.post(f"Removing existing directory {outfile}", "WARN")
                shutil.rmtree(outfile)
            else:
                casalog.post(f"Removing existing file {outfile}", "WARN")
                os.remove(outfile)
        else:
            raise RuntimeError(
                "Logic Error: shouldn't have gotten to this point with overwrite=False"
            )
    with open(outfile, "w") as f:
        json.dump(antpos, f)
