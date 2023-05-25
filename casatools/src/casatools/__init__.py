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
from casaconfig import pull_data, measures_update, config

sakura( ).initialize_sakura( )    ## sakura requires explicit initialization


user_datapath = config.datapath
user_nogui = config.nogui
user_agg = config.agg
user_pipeline = config.pipeline
user_cachedir = __os.path.abspath(__os.path.expanduser(config.cachedir))

user_measurespath = config.measurespath
if (__sys.argv[0] != '-m') and (not __os.path.isdir(__os.path.realpath(user_measurespath))):
    print('measurespath path not found, creating %s' % user_measurespath)
    __os.makedirs(user_measurespath)

ctsys = __utils( )
ctsys.initialize( __sys.executable, user_measurespath, user_datapath, user_nogui,
                  user_agg, user_pipeline, user_cachedir )

logger = logsink(config.logfile) if (hasattr(config,'logfile') and config.logfile is not None) else None

if (__sys.argv[0] != '-m') and (hasattr(config,'configrc')):
    print('Using %s' % config.configrc)
    if logger is not None:
        logger.post('Using %s' % config.configrc, 'INFO')

if hasattr(config,'measures_update') and (config.measures_update == True) and (__sys.argv[0] != '-m'):
    measures_update(ctsys.rundata(), logger=logger)

if __sys.argv[0] != '-m':
    __resolved_iers = ctsys.resolve('geodetic/IERSeop2000')
    if __resolved_iers == 'geodetic/IERSeop2000':
        raise ImportError('measures data is not available, visit https://casadocs.readthedocs.io/en/stable/notebooks/external-data.html for more information')
    if len(ctsys.rundata( )) == 0:
        ctsys.setrundata(__resolved_iers[:-21])

from .coercetype import coerce as __coerce

__coerce.set_ctsys(ctsys)         ## used to locate files from a partial path

def version( ): return list(ctsys.toolversion( ))
def version_string( ): return ctsys.toolversion_string( )

import atexit as __atexit
__atexit.register(ctsys.shutdown) ## c++ shutdown
