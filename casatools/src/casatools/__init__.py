from __future__ import absolute_import
__name__ = 'casatools'
__all__ = [ "ctsys", "version", "version_string"
            'image',
            'logsink',
            'coordsys',
            'synthesisutils',
            'synthesisnormalizer',
            'calanalysis',
            'mstransformer',
            'calibrater',
            'functional',
            'table',
            'measures',
            'imagepol',
            'simulator',
            'sdm',
            'synthesisimstore',
            'miriadfiller',
            'ms',
            'vpmanager',
            'synthesisdeconvolver',
            'vlafiller',
            'sakura',
            'linearmosaic',
            'tablerow',
            'iterbotsink',
            'sidebandseparator',
            'imagemetadata',
            'atcafiller',
            'agentflagger',
            'synthesismaskhandler',
            'regionmanager',
            'msmetadata',
            'imager',
            'singledishms',
            'atmosphere',
            'quanta',
            'synthesisimager',
            'componentlist',
            'spectralline',
            'wvr'
          ]
from .image import image
from .logsink import logsink
from .coordsys import coordsys
from .synthesisutils import synthesisutils
from .synthesisnormalizer import synthesisnormalizer
from .calanalysis import calanalysis
from .mstransformer import mstransformer
from .calibrater import calibrater
from .functional import functional
from .table import table
from .measures import measures
from .imagepol import imagepol
from .simulator import simulator
from .sdm import sdm
from .synthesisimstore import synthesisimstore
from .miriadfiller import miriadfiller
from .ms import ms
from .vpmanager import vpmanager
from .wvr import wvr
from .synthesisdeconvolver import synthesisdeconvolver
from .vlafiller import vlafiller
from .sakura import sakura
from .linearmosaic import linearmosaic
from .tablerow import tablerow
from .iterbotsink import iterbotsink
from .sidebandseparator import sidebandseparator
from .imagemetadata import imagemetadata
from .atcafiller import atcafiller
from .agentflagger import agentflagger
from .synthesismaskhandler import synthesismaskhandler
from .regionmanager import regionmanager
from .msmetadata import msmetadata
from .imager import imager
from .singledishms import singledishms
from .atmosphere import atmosphere
from .quanta import quanta
from .synthesisimager import synthesisimager
from .componentlist import componentlist
from .spectralline import spectralline
from .utils import utils as __utils
import os as __os
import sys as __sys
from casaconfig import get_data_info, do_auto_updates, config
from casaconfig import UnsetMeasurespath, AutoUpdatesNotAllowed, BadLock, BadReadme, NoReadme, RemoteError, NoNetwork
# useful to use here
from casaconfig.private.print_log_messages import print_log_messages

sakura( ).initialize_sakura( )    ## sakura requires explicit initialization


user_datapath = config.datapath
user_nogui = config.nogui
user_agg = config.agg
user_pipeline = config.pipeline
user_cachedir = __os.path.abspath(__os.path.expanduser(config.cachedir))
user_measurespath = config.measurespath
user_verbose = config.casaconfig_verbose

logger = logsink(config.logfile) if (hasattr(config,'logfile') and config.logfile is not None) else None

# data checks, only if user_measurespath is not None
data_info = None
measures_found = False
isSevere = False
# accumuale messages, less confusing when the logger is being redirected to the terminal
msgs = ['']

config_except = None

# first, attempt any auto updates
try:
    # this uses config.measurespath, config.measures_auto_update and config.data_auto_update as appropriate
    do_auto_updates(config, logger)
    
except UnsetMeasurespath as exc:
    msgs.append(str(exc))
    isSevere = True
    msgs.append('')
    msgs.append('Either set measurespath in your personal config.py in ~/.casa or a site config file.')
    msgs.append('CASA may still work if the IERS data can be found in datapath, but this problem is likely to cause casatools to fail to load.')
    msgs.append('')
    # ctsys initialize needs a string for measurespath, leave it empty, it might still work (probably not)
    user_measurespath = ""
    # make sure this prints
    user_verbose = 2
    config_except = exc

except AutoUpdatesNotAllowed as exc:
    msgs.append(str(exc))
    msgs.append('')
    msgs.append('Warning: measurespath must exist as a directory and it must be owned by the user.')
    msgs.append('Warning: no auto update is possible on this measurespath by this user.')
    msgs.append('')
    # this is reraised only if the data isn't found in datapath
    # make sure this prints
    user_verbose = 2
    config_except = exc

except BadLock as exc:
    msgs.append(str(exc))
    # this possibly indicates a serious problem, reraise this if the data can't be found
    # it's severe, it will print
    isSevere = True
    config_except = exc

except BadReadme as exc:
    msgs.append(str(exc))
    # this likely indicates a problem, reraise this if the data can't be found
    msgs.append('')
    msgs.append('Some or all of the expected auto updates did not happen.')
    msgs.append('This indicates something went wrong on a previous update and the data should be reinstalled.')
    msgs.append('Updates will continue to fail until the data are reinstalled.')
    msgs.append('If the IERSeop2000 table is found in datapath then casatools will import.')
    msgs.append('')
    # it's severe, it will print
    isSevere = True
    config_except = exc

except NoReadme as exc:
    # this is a symptom that's fully explained here and possibly again if this is reraised, don't print out str(exc) here
    msgs.append('Some or all of the expected auto updates did not happen.')
    msgs.append('This indicates that measurespath is not empty and does not contain data maintained by casaconfig.')
    msgs.append('If the IERSeop2000 table is found in datapth then casatools will import.')
    msgs.append('')
    # print this
    user_verbose = 2
    config_except = exc

except NoNetwork as exc:
    msgs.append('')
    msgs.append('No data or measures updates could be done because there is no network connection.')
    # print it, this is ONLY reraised if the IERSeop2000 table is not found so it's clearer what the root cause likely was.
    user_verbose = 2
    config_except = exc

except RemoteError as exc:
    msgs.append(str(exc))
    msgs.append('')
    msgs.append('Some or all of the expepcted auto updates did not happen.')
    msgs.append('Either there is no network connection, there is no route to the remote server, or the remote server is offline')
    msgs.append('If the IERSeop2000 table is found in datapath then casatools will import without any updates. Try again later for updates')
    msgs.append('')
    # print this
    user_verbose = 2
    config_except = exc
    
except Exception as exc:
    msgs.append('ERROR! Unexpected exception while doing auto updates or checking on the status of the data at measurespath')
    print(str(exc))
    # print this
    # this seems severe, it will print
    isSevere = True
    config_except = exc

# now try and get any info
# skip this step if UnsetMeasurespath was previously raised
if not isinstance(config_except,UnsetMeasurespath):
    try:
    
        data_info = get_data_info(user_measurespath, logger)
        if data_info['casarundata'] is None:
            isSevere = True
            msgs.append('The expected casa data was not found at measurespath. CASA may still work if the data can be found in datapath.')
        elif data_info['casarundata']['version'] == 'invalid':
            isSevere = True
            msgs.append('The contents of measurespath do not appear to be casarundata. CASA may still work if the data can be found in datapath.')
        elif data_info['casarundata']['version'] == 'unknown':
            msgs.append('The casa data found at measurespath is not being maintained by casaconfig.')
        else:
            pass

        if data_info['measures'] is None:
            isSevere = True
            msgs.append('The expected measures data was not found at measurespath. CASA may still work if the data can be found in datapath.')
        elif data_info['measures']['version'] == 'invalid':
            isSevere = True
            msgs.append('The contents of measurespath do not appear to include measures data. CASA may still work if the data can be found in datapath.')
        elif data_info['measures']['version'] == 'unknown':
            msgs.append('The measures data found at measurespath is not being maintained by casaconfig.')
            measures_found = True
        else:
            measures_found = True

    except Exception as exc:
        # no exceptions are expected here
        msgs.append('ERROR! Unexpected exception while checking on the status of the data at measurespath')
        # this seems severe
        isSevere = True
        print(str(exc))
        config_except = exc

if (len(msgs) > 1):
    msgs.append('visit https://casadocs.readthedocs.io/en/stable/notebooks/external-data.html for more information')
    msgs.append('')
    print_log_messages(msgs, logger, isSevere, verbose=user_verbose)

# don't use user_measurespath here if not ok
if not measures_found:
    user_measurespath = ""
else:
    # always append measurespath to datapath if not already found there to be used by ctsys.initialize
    add_mp = True
    for apath in user_datapath:
        if __os.path.exists(apath) and __os.path.samefile(apath, user_measurespath):
            add_mp = False
            break
    if add_mp:
        user_datapath.append(user_measurespath)

ctsys = __utils( )
ctsys.initialize( __sys.executable, user_measurespath, user_datapath, user_nogui,
                  user_agg, user_pipeline, user_cachedir )

# try and find the IERS data
__resolved_iers = ctsys.resolve('geodetic/IERSeop2000')
if __resolved_iers == 'geodetic/IERSeop2000':
    # this match means ctsys.resolve did not found it in datapath
    # if there was a previously raised casaconfig exception, re-raise it here so that casashell can know what went wrong
    if config_except is not None:
        print("measures data is not available")
        raise config_except
    else:
        raise ImportError('measures data is not available, visit https://casadocs.readthedocs.io/en/stable/notebooks/external-data.html for more information')
else:
    # if measures was not previously found, use this path as measurespath
    if not measures_found:
        # this removes the "geodetic/IERSeop2000" and the preceding "/" from the returned path, that's the measurespath to be used here
        
        user_measurespath = __resolved_iers[:-21]
        ctsys.setmeasurespath(user_measurespath)
        print_log_messages(["Using location of IERSeop2000 table found in datapath for measurespath %s" % user_measurespath], logger)
    else:
        # check that this path is the same as what is expected at measurespath
        if not __os.path.samefile(__resolved_iers, __os.path.join(user_measurespath,'geodetic/IERSeop2000')):
            print("\nThe path to the geodetic IERSeop2000 table in measurespath is not the same location found first in datapath for that table.")
            print("CASA should work in this configuration. The datapath list will be used to search for data needed by CASA.")
            print("\nThe measurespath value is used to find that IERS table and this indicates that the")
            print("measures tables present in datapath may be different from those found in measurespath.")
            print("\nIf this was not expected you may want to use your config file (normally at ~/.casa/config.py) to")
            print("set datapath to put measurespath first or set measurespath to include the measures data found in datapath.")
            logger.post("WARNING: geodetic/IERSeop2000 found at measurespath is not the same table as found in datapath",'WARN')

from .coercetype import coerce as __coerce

__coerce.set_ctsys(ctsys)         ## used to locate files from a partial path

def version( ): return list(ctsys.version( ))
def version_string( ): return ctsys.version_string( )

import atexit as __atexit
__atexit.register(ctsys.shutdown) ## c++ shutdown
