import os

from casatools import ms
_ms = ms( )

def listfits(fitsfile=None):
    """

    """
    #Python script

    if ((type(fitsfile)==str) & (os.path.exists(fitsfile))):
        _ms.listfits(fitsfile);
    else:
        raise ValueError('fits file not found - please verify the name')
