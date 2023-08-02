//# WPConvFunc.h: Definition for WPConvFunc
//# Copyright (C) 2007-2016
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be adressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#
//# $Id$

#ifndef SYNTHESIS_TRANSFORM2_WPCONVFUNC_H
#define SYNTHESIS_TRANSFORM2_WPCONVFUNC_H

#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Containers/Block.h>
#include <casacore/casa/Utilities/CountedPtr.h>
#include <casacore/casa/Arrays/ArrayFwd.h>

namespace casacore{

	template<class T> class ImageInterface;

}

namespace casa{
	namespace vi{class VisBuffer2;}

namespace refim{// namespace for imaging refactoring

  // <summary>  A class to support FTMachines get their convolution casacore::Function </summary>
  
  // <use visibility=export>
  // <prerequisite>
  //   <li> <linkto class=VisBuffer>VisBuffer</linkto> module
// </prerequisite>
  // <etymology>
  // WP for W-casacore::Projection 
  // ConvFunc => returns the convolution functions
  // </etymology>
  //
  // <synopsis> 
  // FTMachines like WProjection and MosaicFT need convolution functions to 
  // deal with directional dependent issues...
  // this class and related ones provide and cache  such functions for re-use 
  //</synopsis>

  class WPConvFunc 
    {
    public:
      WPConvFunc(const casacore::Double minW=-1.0, const casacore::Double maxW=-1.0, const casacore::Double rmsW=-1.0);
      WPConvFunc(const casacore::RecordInterface& rec);
      //Copy constructor
      WPConvFunc(const WPConvFunc& other);
      //
      WPConvFunc& operator=(const WPConvFunc&other);
      

      virtual ~WPConvFunc();

      // Inputs are the image, visbuffer,  wConvsize
      // findconv return a cached convolution function appropriate for this 
      // visbuffer and number of w conv plane
      void findConvFunction(const casacore::ImageInterface<casacore::Complex>& iimage, 
			    const vi::VisBuffer2& vb,
			    const casacore::Int& wConvSize,
			    const casacore::Vector<casacore::Double>& uvScale,
			    const casacore::Vector<casacore::Double>& uvOffset,
			    const casacore::Float& padding, 
			    casacore::Int& convSampling,
			    casacore::Cube<casacore::Complex>& convFunc, 
			    casacore::Int& convsize,
			    casacore::Vector<casacore::Int>& convSupport,
			    casacore::Double& wScale);

      
      virtual casacore::Bool makeAverageResponse(const vi::VisBuffer2& /*vb*/,
				       const casacore::ImageInterface<casacore::Complex>& /*image*/,
				     //				     casacore::TempImage<casacore::Float>& theavgPB,
				       casacore::ImageInterface<casacore::Float>& /*theavgPB*/,
				       casacore::Bool /*reset=true*/)
    {throw(casacore::AipsError("WPConvFunc::makeAverageRes() called"));};
    ///Make full WConfFunction, despite it being circularly symmetric; can be used along 
    // with A-term convolution for a Vector of W values
    // is the coordinateSystem to get the scale of pixels
    // csys is the image based csys it will be returned in the UV domain
    casacore::Bool makeWConvFuncs(casacore::Cube<casacore::Complex>& wconv, casacore::Vector<casacore::Int>& supports,  casacore::CoordinateSystem& cs, const casacore::Int& npix, const casacore::Vector<casacore::Double>& wVals); 
    
    // wVal is the w-value in lambda
    casacore::Bool makeSkyWFunc(casacore::Matrix<casacore::Complex>& wSkyFunc, const casacore::CoordinateSystem& cs, const casacore::Int& npix, const casacore::Double& wVal); 
    
      //Serialization
      casacore::Bool toRecord(casacore::RecordInterface& rec);
      casacore::Bool fromRecord(casacore::String& err, const casacore::RecordInterface& rec);
    private:
      casacore::Bool checkCenterPix(const casacore::ImageInterface<casacore::Complex>& image);
      void makeGWplane(casacore::Matrix<casacore::Complex>& screen, const casacore::Int iw, casacore::Double s0, casacore::Double s1, casacore::Float *& wsaveptr, casacore::Int& lsav, casacore::Int& inner, casacore::Complex*& cor, casacore::Double&cpWscale);
      casacore::Int findSupport(casacore::Matrix<casacore::Complex>& scr); 
      casacore::Block <casacore::CountedPtr<casacore::Cube<casacore::Complex> > > convFunctions_p;
      casacore::Block <casacore::CountedPtr<casacore::Vector<casacore::Int> > > convSupportBlock_p;
      std::map <casacore::String, casacore::Int> convFunctionMap_p;
      casacore::Vector<casacore::Int> convSizes_p;

      casacore::Int actualConvIndex_p;
      casacore::Int convSize_p;
      casacore::Vector<casacore::Int> convSupport_p;
      casacore::Cube<casacore::Complex> convFunc_p;
      casacore::Double wScaler_p;
      casacore::Int convSampling_p;
      casacore::Int nx_p, ny_p;
      casacore::Double minW_p, maxW_p, rmsW_p;

    };
} //end of namespace refim
} // end namespace casa
#endif
