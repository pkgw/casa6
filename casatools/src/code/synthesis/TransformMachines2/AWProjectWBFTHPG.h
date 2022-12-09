//# AWProjectFT.h: Definition for AWProjectFT
//# Copyright (C) 1996,1997,1998,1999,2000,2002
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
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

#ifndef SYNTHESIS_TRANSFORM2_AWPROJECTWBFTHPG_H
#define SYNTHESIS_TRANSFORM2_AWPROJECTWBFTHPG_H
#include <synthesis/TransformMachines2/FTMachine.h>
#include <synthesis/TransformMachines2/AWProjectWBFT.h>

namespace casa
{ //# NAMESPACE CASA - BEGIN
  namespace refim
  {
    class AWProjectWBFTHPG : public AWProjectWBFT
    {
    public:
      AWProjectWBFTHPG(casacore::Int nFacets, casacore::Long cachesize,
		       casacore::CountedPtr<CFCache>& cfcache,
		       casacore::CountedPtr<ConvolutionFunction>& cf,
		       casacore::CountedPtr<VisibilityResamplerBase>& visResampler,
		       casacore::Bool applyPointingOffset=true,
		       vector<float> pointingOffsetSigDev = {10,10},
		       casacore::Bool doPBCorr=true,
		       casacore::Int tilesize=16, 
		       casacore::Float paSteps=5.0, 
		       casacore::Float pbLimit=5e-4,
		       casacore::Bool usezero=false,
		       casacore::Bool conjBeams_p=true,
		       casacore::Bool doublePrecGrid=false):
	AWProjectWBFT(nFacets, cachesize, cfcache, cf, visResampler,
		      applyPointingOffset, pointingOffsetSigDev,
		      doPBCorr, tilesize, paSteps, pbLimit, usezero,
		      conjBeams_p, doublePrecGrid),applyFFT_p(false)
      {};
      
      ~AWProjectWBFTHPG(){};
      
      // Assignment operator
      AWProjectWBFTHPG &operator=(const AWProjectWBFTHPG &other)
      {
	if(this!=&other) 
	  {
	    //Do the base parameters
	    AWProjectWBFT::operator=(other);
	    
	    applyFFT_p=other.applyFFT_p;
	  }
	return *this;
      };
      
      //---------------------------------------------------------------------------------------
      // Overloading getImage() to not do FFT here.  The appropriate
      // FFT is applied on the GPU.  Here, only conversion and copy
      // from DP to SP image is done.
      //
      virtual casacore::ImageInterface<casacore::Complex>&
      getImage(casacore::Matrix<casacore::Float>& weights,
	       casacore::Bool normalize=false);

      virtual void getWeightImage(casacore::ImageInterface<casacore::Float>& weightImage, casacore::Matrix<casacore::Float>& weights);

      virtual casacore::String name() const { return "AWProjectWBFTHPG";};

      virtual void resampleDataToGrid(casacore::Array<casacore::Complex>& griddedData,VBStore& vbs,
				      const VisBuffer2& vb, casacore::Bool& dopsf);
      virtual void resampleDataToGrid(casacore::Array<casacore::DComplex>& griddedData,VBStore& vbs,
				      const VisBuffer2& vb, casacore::Bool& dopsf);


      ///re implement the initializetoVis as there is no FFT needed for model on gpu
      virtual void initializeToVisNew(const vi::VisBuffer2& vb,
					     casacore::CountedPtr<SIImageStore> imstore);
    private:
      
      Bool applyFFT_p;
    };
  } //# NAMESPACE CASA - END
};
#endif
