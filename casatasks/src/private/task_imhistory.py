
from casatools import image
from .. import casalog

def imhistory(
    imagename, mode, verbose, origin, message
):
    _myia = image()
    try:
        casalog.origin('imhistory')
        _myia.open(imagename)
        if mode.startswith("l") or mode.startswith("L"):
            return _myia.history(verbose)
        elif mode.startswith("a") or mode.startswith("A"):
            return _myia.sethistory(origin=origin, history=message)
        raise ValueError("Unsopported mode " + mode)

    finally:
        _myia.done()
