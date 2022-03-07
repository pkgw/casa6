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
#*  Created on: Mar 07, 2022
#*      Author: kgolap
#*
from __future__ import absolute_import

import os
import shutil
# get is_CASA6 and is_python3
from casatasks.private.casa_transition import *

from casatasks import casalog
from casatools import msuvbinner as msbin
from casatools import msmetadata
from casatools import ms
ms=ms()
msmd=msmetadata()
from  casatasks.private.imagerhelpers.input_parameters import saveparams2last

@saveparams2last(multibackup=True) 
def msuvbin(vis=None, field=None, spw=None, taql=None, outvis=None, phasecenter=None, nx=None, ny=None, cell=None,
            ncorr=None, nchan=None, fstart=None, fstep=None, wproject=None, memfrac=None, doflag=None):
    
    casalog.origin('msuvbin ')
    if(field==''):
        field='*'
    fieldid=0
    fieldid=ms.msseltoindex(vis=vis, field=field)['field'][0]
    if(phasecenter==''):
        msmd.open(vis)
        phcen=msmd.phasecenter(fieldid)
        msmd.done()
        phasecenter=phcen['refer']+' '+str(phcen['m0']['value'])+str(phcen['m0']['unit'])+' '+str(phcen['m1']['value'])+str(phcen['m1']['unit'])
    msbinner=msbin(phasecenter=phasecenter, nx=nx, ny=ny, ncorr=ncorr, nchan=nchan, cellx=cell, celly=cell, fstart=fstart, fstep=fstep, memfrac=memfrac, wproject=wproject, doflag=doflag)
    msbinner.selectdata(msname=vis, spw=spw, field=field, taql=taql)
    msbinner.setoutputms(outvis)
    msbinner.filloutputms()
