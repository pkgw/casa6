"""
 This is the private module within casatasks which mediates access to
 user configuration options. It selects which user configuration file
 to load (only one should be loaded; if desired that file should
 load secondary files). config.py then pulls supported flags from from
 this file after validating them.
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
exec(open(configrc).read())

# call casashell's casataskrc in case there is a command line argument to overwrite the rcdir
try:
    from casataskrc import *
except:
    pass

