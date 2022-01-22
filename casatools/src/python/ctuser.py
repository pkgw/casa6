"""Hook to allow user-specified customization code to run.

This code was imported from Python 2.7 for use with CASAtools...

However, some programs or sites may find it convenient to allow users
to have a standard customization file, which gets run when a program
requests it.  This module implements such a mechanism.  A program
that wishes to use the mechanism must execute the statement

    import ctuser

"""
import os
import pkg_resources

home = os.curdir  # Default
if 'HOME' in os.environ:
    home = os.environ['HOME']
elif os.name == 'posix':
    home = os.path.expanduser("~/")
elif os.name == 'nt':  # Contributed by Jeff Bauer
    if 'HOMEPATH' in os.environ:
        if 'HOMEDRIVE' in os.environ:
            home = os.environ['HOMEDRIVE'] + os.environ['HOMEPATH']
        else:
            home = os.environ['HOMEPATH']

# use ~/.casa/config.py if it exists
# otherwise use default from casaconfig package
configrc = os.path.join(home, ".casa/config.py")
configrc = configrc if os.path.exists(configrc) else pkg_resources.resource_filename('casaconfig', 'config.py')

# call casashell's casatoolrc in case there is a command line argument to overwrite the rcdir
try:
    exec(open(configrc).read())
    from casatoolrc import *
except:
    pass
