import os

from casatools import table
_tb = table( )

def clearstat():
    """Clear all read/write locks on tables. This can be used if a task has
       indicated that it is trying to get a lock on a file.

    """
    _tb.clearlocks( )
