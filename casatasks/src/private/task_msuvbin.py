# msuvbin task
# Copyright (C) 2022
# Associated Universities, Inc. Washington DC, USA.
#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
# License for more details.
#
# https://www.gnu.org/licenses/
#
# Queries concerning CASA should be submitted at
#        https://help.nrao.edu
#
#        Postal address: CASA Project Manager
#                        National Radio Astronomy Observatory
#                        520 Edgemont Road
#                        Charlottesville, VA 22903-2475 USA
#
# $Id$
# *  Created on: Mar 07, 2022
# *      Author: kgolap
# *

import os
import shutil
import typing
from typing import Tuple, List, Union, Optional

from . import flaghelper as fh

from casatasks import casalog
from casatools import msuvbinner as msbin
from casatools import msmetadata
from casatools import ms

ms = ms()
msmd = msmetadata()
from casatasks.private.imagerhelpers.input_parameters import (
    saveparams2last,
    determineFreqRange,
)


@saveparams2last(multibackup=True)
def msuvbin(
    vis: Optional[str] = None,
    field: Optional[str] = None,
    spw: Optional[str] = None,
    taql: Optional[str] = None,
    outputvis: Optional[str] = None,
    phasecenter: Optional[str] = None,
    imsize: Optional[Union[List[int], List[float], int, float]] = None,
    cell: Optional[str] = None,
    ncorr: Optional[int] = None,
    nchan: Optional[int] = None,
    start: Optional[str] = None,
    width: Optional[str] = None,
    wproject: Optional[bool] = None,
    memfrac: Optional[float] = None,
    mode: Optional[str] = None,
    flagbackup: Optional[bool] = None,
) -> None:
    fstart = start
    fstep = width
    casalog.origin("msuvbin ")
    if wproject:
        casalog.post(
            "The wprojection option is extremely slow; you may consider running without it",
            "WARN",
            "task_msuvbin"
        )
    if field == "":
        field = "*"
    fieldid = 0
    fieldid = ms.msseltoindex(vis=vis, field=field)["field"][0]
    if isinstance(imsize, (int, float)):
        nx = imsize
        ny = imsize
    else:
        nx = imsize[0]
        ny = imsize[0] if (len(imsize)==1) else imsize[1]
    if phasecenter == "":
        msmd.open(vis)
        phcen = msmd.phasecenter(fieldid)
        msmd.done()
        phasecenter = (
            phcen["refer"]
            + " "
            + str(phcen["m0"]["value"])
            + str(phcen["m0"]["unit"])
            + " "
            + str(phcen["m1"]["value"])
            + str(phcen["m1"]["unit"])
        )
    if spw == "":
        spw = "*"
    if(nchan < 1):
        casalog.post(
            "nchan has to be larger than 0", "ERROR",
            "task_msuvbin"
        )
        
    if (not start) or (start == ""):
        (fbeg, fwidth) = determineFreqRange(vis=vis, fieldid=fieldid, spw=spw)
        fstart = f"{fbeg}Hz"
        if (not width) or (width == ""):
            fstep = f"{fwidth/nchan}Hz"
    #print(f"fstart={fstart}, fstep={fstep}")
    doflag= "write_flags" in mode 
    if(doflag and flagbackup):
        fh.backupFlags(aflocal=None, msfile=vis, prename='msuvbin')
    msbinner = msbin.msuvbinner(
        phasecenter=phasecenter,
        nx=nx,
        ny=ny,
        ncorr=ncorr,
        nchan=nchan,
        cellx=cell,
        celly=cell,
        fstart=fstart,
        fstep=fstep,
        memfrac=memfrac,
        wproject=wproject,
        doflag=doflag,
    )
    msbinner.selectdata(msname=vis, spw=spw, field=field, taql=taql)
    msbinner.setoutputms(outputvis)
    msbinner.filloutputms()
    del msbinner
