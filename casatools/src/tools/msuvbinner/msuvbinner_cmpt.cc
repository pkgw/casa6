/*
//# msuvbinner_cmpt tool for binding to python 
//# Copyright (C) 2022
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by
//# the Free Software Foundation; either version 3 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
//# License for more details.
//#
//# https://www.gnu.org/licenses/
//#
//# Queries concerning CASA should be submitted at
//#        https://help.nrao.edu
//#
//#        Postal address: CASA Project Manager 
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$
 *  Created on: Feb 28, 2022
 *      Author: kgolap
 */


#include <casacore/casa/Quanta/QuantumHolder.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>
#include <casacore/coordinates/Coordinates/Projection.h>
#include <casacore/measures/Measures/MDirection.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/VisibilityIterator2.h>
#include <casacore/ms/MeasurementSets/MSMainColumns.h>
#include <mstransform/MSTransform/MSTransformDataHandler.h>
#include <mstransform/MSTransform/MSUVBin.h>
#include <casacore/tables/Tables/TableIter.h>
#include <casacore/tables/Tables/TableColumn.h>
#include <casacore/casa/Inputs/Input.h>
#include <casacore/casa/namespace.h>
#include <msuvbinner_cmpt.h>
using namespace std;
using namespace casa;
using namespace casacore;
namespace casac {

	
	msuvbinner::msuvbinner(const ::casac::variant& phasecenter, long nx, long ny, long ncorr, long nchan, const string& cellx, const string& celly, const string& fstart, const string& fstep, double memfrac, bool wproject, bool doflag ){
		 itsLog = new LogIO();
		 try{
			casacore::MDirection phaseCenter;
			if(!casaMDirection(phasecenter, phaseCenter)){
				*itsLog << "Could not interprete phasecenter param : "+phasecenter.toString() << LogIO::SEVERE <<LogIO::POST;
			}
			casacore::Quantity cellX=casaQuantity(cellx);
			casacore::Quantity cellY=casaQuantity(celly);
			casacore::Quantity freqStart=casaQuantity(fstart);
			casacore::Quantity freqStep=casaQuantity(fstep);
  
			itsBinner=MSUVBin(phaseCenter, nx,
				ny, nchan, ncorr, cellX, cellY, freqStart, freqStep, memfrac, wproject, doflag);
			 
		 }
		 catch(AipsError x){
			RETHROW(x); 
		 }
		
	}
	msuvbinner::~msuvbinner(){
		
	}
	////////////////////////////
	bool msuvbinner::selectdata(const string& msname, const string& spw, const string& field, const string& baseline, const string& scan, const string& uvrange, const string& taql){
		try{
			itsBinner.selectData(msname, spw, field, baseline, scan, uvrange, taql);
			return true;
		}
		catch(AipsError x){
			RETHROW(x); 
		}
		 
		return false;
		
	}
	///////////////////////////////////
	bool msuvbinner::setoutputms(const string& outms){
		try{
			itsBinner.setOutputMS(outms);
			return true;
		}
		catch(AipsError x){
			RETHROW(x); 
		}
		 
		return false;	
		
		
	}
	
	bool msuvbinner::filloutputms(){
		try{
			return itsBinner.fillOutputMS();

		}
		catch(AipsError x){
			RETHROW(x); 
		}
		 
		return false;	
		
		
	}
	
	
}//namespace casac
	
	
