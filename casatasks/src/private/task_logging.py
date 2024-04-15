from casatasks import casalog as _clog
from datetime import datetime as _time
import casatasks
import os 

def start_log( tname, arguments ):
    spaces = ' '*(18-len(tname))
    start_time = str(_time.now())
    _clog.origin(tname)
    _clog.post( '##########################################' )
    _clog.post( '##### Begin Task: ' + tname + spaces + ' #####' )
    _clog.post( '%s( %s )' % ( tname, ', '.join(arguments) ))
    return start_time,

def end_log( state, tname, result ):
    spaces = ' '*(18-len(tname))
    end_time = str(_time.now())
    _clog.origin(tname)
    _clog.post( 'Result {}: {}'.format(tname, repr(result)), priority='DEBUG')
    _clog.post( 'Task ' + tname + ' complete. Start time: ' + state[0] + ' End time: ' + end_time )
    _clog.post( '##### End Task: ' + tname + '  ' + spaces + ' #####' )
    _clog.post( '##########################################' )
    return result

def except_log(tname, exc):
    _clog.post('Task {} raised an exception of class {} with the following message: {}'.
               format(tname, type(exc).__name__, exc) ,'ERROR')
