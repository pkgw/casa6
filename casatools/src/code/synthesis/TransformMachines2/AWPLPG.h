//# AWPLPG.h: Definition for a CPU based gridder for A,W (LPG= LowPerformanceGridder)
//# Copyright (C) 2023
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
//# You should have received a copy of the GNU  General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Queries concerning CASA should be submitted at
//#        https://help.nrao.edu
//#
//#        Postal address: CASA Project Manager 
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#

#ifndef SYNTHESIS_TRANSFORM2_AWPLPG_H
#define SYNTHESIS_TRANSFORM2_AWPLPG_H

#include <synthesis/TransformMachines2/MosaicFTNew.h>

namespace casa { //# NAMESPACE CASA - BEGIN

namespace refim{
  class AWConvFuncHolder;
/*
  An copy MosaicFT except
  Looks like it is just to get  differentlently normalized images i.e (image*nx*ny)
  which implies somewhere  some code  is   just be using FFT in forward (the toFrequency direction) 
  when it should be the reverse (or vice-versa)
*/

class AWPLPG : public MosaicFTNew {
public:

 

  AWPLPG(SkyJones* sj, const casacore::Int nw,  const casacore::Bool dosquint, const casacore::Double painc, casacore::MPosition mloc, casacore::String stokes,  const casacore::Bool usezero=true, const casacore::Bool useDoublePrec=true,  const casacore::Bool usePointing=false);
  AWPLPG(const AWPLPG &other);
  AWPLPG& operator=(const AWPLPG& other);
  virtual refim::FTMachine* cloneFTM();
    // Get actual coherence from grid by degridding
  virtual void get(vi::VisBuffer2& vb, casacore::Int row=-1);


  // Put coherence to grid by gridding.
  virtual void put(const vi::VisBuffer2& vb, casacore::Int row=-1, casacore::Bool dopsf=false, 
	   FTMachine::Type type=FTMachine::OBSERVED);

  virtual void gridImgWeights(const vi::VisBuffer2& vb);
  
  
protected:     
  
  virtual void findConvFunction(const casacore::ImageInterface<casacore::Complex>& image,
			const vi::VisBuffer2& vb);
  virtual void init(const vi::VisBuffer2& vb);
  
  std::shared_ptr<AWConvFuncHolder> awConvs_p;
  casacore::Bool doSquint_p;
  casacore::Double paInc_p;
  casacore::Int nw_p;
  
};
  } //refim ends
} //# NAMESPACE CASA - END

#endif
