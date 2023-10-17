
//# AWConvFuncHolder.h: Definition for a holder class for conv functions for mosaic with squint and wproject
//# Copyright (C) 2018-2019
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

#ifndef SYNTHESIS_TRANSFORM2_AWCONVFUNCHOLDER_H
#define SYNTHESIS_TRANSFORM2_AWCONVFUNCHOLDER_H
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <synthesis/TransformMachines2/AWConvFunc.h>
#include <msvis/MSVis/VisBufferUtil.h>
namespace casa{ //# namespace casa
namespace refim{ //#	 namespace for refactored imaging code with vi2/vb2
class AWConvFuncHolder{
  public:
  AWConvFuncHolder(const casacore::CoordinateSystem& csys, const int nx, const int ny, 
                   const bool dosquint=false, const double paInc=2*C::pi, const casacore::String& obs="EVLA", const int& oversamp=4);
  AWConvFuncHolder(const AWConvFuncHolder& other);
  AWConvFuncHolder& operator=(const AWConvFuncHolder& other);
  bool loadConvFromFile();
  bool saveConvToFile(const casacore::String& filename);
  bool addConvFunc(const casacore::Vector<casacore::Double>& freqs, const casacore::Vector<Double>& wVals,  const casacore::Double& painc);
  //Get the convFunctions and indexing from a given vb
  casacore::Array<casacore::Complex>& getConvFunc();
  casacore::Array<casacore::Complex>& getWeightConvFunc();
  casacore::Vector<casacore::Int> getConvSizes();
  casacore::Vector<casacore::Int> getConvSupports();
  casacore::Vector<Stokes::StokesTypes> getPolVals(){return polVals_p;};
  casacore::Vector<Double> getFreqVals(){return freqVals_p;};
  casacore::Vector<Double> getWVals(){return wVals_p;};
  casacore::Vector<Double> getPAVals(){return paVals_p;};
  int getOverSampling(){return oversamp_p;};
  //Rowmap will return the indices to match along the 5th axis of convFunc, polmap is for the 3rd axis, and chanmap is for the 4th axis.
  //Rowmap will map combination of pa, antennapair and w to give the 5th index that matches 
  //Rowmap will be the same nrow as vb.nrows , polmap will gave the same length of vb.ncorrelations and chanmap will be the length of vb.nchannelscasacore::Vector<casacore::Int>& rowMap
  void getConvIndices( casacore::Vector<casacore::Int>& polMap, casacore::Vector<casacore::Int>& chanMap, casacore::Vector<casacore::Int>& rowMap, const vi::VisBuffer2& vb, const casacore::Matrix<casacore::Double>& rotuvw);
  //Function gives pointing direction w.r.t image center  in phase shift: used in putting a phase gradient in the UV domain
  // Will be using for now only 1st row.
  Vector<Double> getPointingPhaseShift(const vi::VisBuffer2& vb, const bool usepointing=False);
  

 private:
   void appendConvFuncs(const casacore::Array<casacore::Complex>& awConv,  const casacore::Array<casacore::Complex>& aWwtConv,  const casacore::Matrix<casacore::Int>& awsupport, const casacore::Vector<casacore::Double>& newfreqs, const casacore::Double paval);
   
   
  double painc_p;
  bool dosquint_p;
  casacore::Array<casacore::Complex> convFunc_p;
  casacore::Array<casacore::Complex> wgtConvFunc_p;
  casacore::Vector<Stokes::StokesTypes> polVals_p;
  casacore::Vector<Double> freqVals_p;
  casacore::Vector<Double> wVals_p;
  casacore::Vector<Double> paVals_p;
  casacore::Vector<std::pair<int, int> > antpairVals_p;
  ///These 3 vectors will be the same length as the 5th axis pf convFunc
  // they will point to the index of wvals or paval or antpair vetcor value for which the convfunc is
  casacore::Vector<Int> rowAxisWVals_p;
  casacore::Vector<Int> rowAxisPAVals_p;
  casacore::Vector<Int > rowAxisAntennaPair_p;
  casacore::Vector<Int> convSizes_p;
  casacore::Vector<Int> convSupport_p;
  //Image parameters that the convfunc will map to (for now the direction
  //increment is what matters for the csys_p
  casacore::CoordinateSystem outcsys_p;
  casacore::CoordinateSystem calcCsys_p;
  int nx_p;
  int ny_p;
  int calcNpix_p;
  int oversamp_p;
  std::shared_ptr<EVLAAperture> aterm_p;
  std::shared_ptr<VisBufferUtil> vbutil_p;                 
  
};
  
   }//# end namespace refim
} // end namespace casa


#endif








