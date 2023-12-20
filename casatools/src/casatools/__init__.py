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
# useful to use here
from casaconfig.private.print_log_messages import print_log_messages

sakura( ).initialize_sakura( )    ## sakura requires explicit initialization


user_datapath = config.datapath
user_nogui = config.nogui
user_agg = config.agg
user_pipeline = config.pipeline
user_cachedir = __os.path.abspath(__os.path.expanduser(config.cachedir))
user_measurespath = config.measurespath

logger = logsink(config.logfile) if (hasattr(config,'logfile') and config.logfile is not None) else None

# this uses config.measurespath, config.measures_auto_update and config.data_auto_update as appropriate
do_auto_updates(config, logger)

# data checks, only if user_measurespath is not None
data_info = None
if user_measurespath is not None:
    data_info = get_data_info(user_measurespath, logger)
    data_ok = False
    if data_info['casarundata'] is None:
        print_log_messages('The expected casa data was not found at measurespath. CASA may still work if the data can be found in datapath.', logger, True)
    elif data_info['casarundata'] == 'invalid':
        print_log_messages('The contents of measurespath do not appear to be casarundata. CASA will likely fail as a result', logger, True)
    elif data_info['casarundata'] == 'unknown':
        print_log_messages('The casa data found at measurespath is not being maintained using casaconfig tools. CASA will still work but that data may be out of date.', logger)
    else:
        data_ok = True

    measures_ok = False
    if data_info['measures'] is None:
        print_log_messages('The expected measures data was not found at measurespath. CASA may still work if the data can be found in datapath.', logger, True)
    elif data_info['casarundata'] == 'invalid':
        print_log_messages('The contents of measurespath do not appear to include measures data. CASA will likely fail as a result', logger, True)
    elif data_info['measures'] == 'unknown':
        print_log_messages('The measures data found at measurespath is not being maintained using casaconfig tools. CASA will still work but that data may be out of date.', logger)
    else:
        measures_ok = True

    if (not data_ok) or (not measures_ok):
        print('visit https://casadocs.readthedocs.io/en/stable/notebooks/external-data.html for more information')

ctsys = __utils( )
ctsys.initialize( __sys.executable, user_measurespath, user_datapath, user_nogui,
                  user_agg, user_pipeline, user_cachedir )

# try and find the IERS data
__resolved_iers = ctsys.resolve('geodetic/IERSeop2000')
if __resolved_iers == 'geodetic/IERSeop2000':
    raise ImportError('measures data is not available, visit https://casadocs.readthedocs.io/en/stable/notebooks/external-data.html for more information')

# and use this as rundata (measurespath) if the data_info for measures is None or the measures version is invalid
if data_info['measures'] is None or data_info['measures']['version'] == "invalid":
    ctsys.setrundata(__resolved_iers[:-21])

from .coercetype import coerce as __coerce

__coerce.set_ctsys(ctsys)         ## used to locate files from a partial path

def version( ): return list(ctsys.toolversion( ))
def version_string( ): return ctsys.toolversion_string( )

import atexit as __atexit
__atexit.register(ctsys.shutdown) ## c++ shutdown
