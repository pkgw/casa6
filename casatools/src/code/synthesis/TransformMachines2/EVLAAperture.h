// -*- C++ -*-
//# EVLAAperture.h: Definition of the EVLAAperture class
//# Copyright (C) 1997,1998,1999,2000,2001,2002,2003
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
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$
//
#ifndef SYNTHESIS_TRANSFORM2_EVLAAPERTURE_H
#define SYNTHESIS_TRANSFORM2_EVLAAPERTURE_H

#include <synthesis/TransformMachines2/Utils.h>
#include <casacore/images/Images/ImageInterface.h>
#include <synthesis/TransformMachines/BeamCalc.h>
//#include <synthesis/MeasurementComponents/ATerm.h>
#include <synthesis/TransformMachines2/AzElAperture.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/casa/Arrays/ArrayFwd.h>
//
//---------------------------------------------------------------------
//

namespace casa { //# NAMESPACE CASA - BEGIN
  
  class VisBuffer;
  //  class EVLAAperture : public ATerm
  namespace refim{
  class EVLAAperture : public AzElAperture
  {
  public:
    //    EVLAAperture(): ATerm(), polMap_p(), feedStokes_p() {};
    EVLAAperture();
    ~EVLAAperture() {};
    EVLAAperture& operator=(const EVLAAperture& other);
    //
    // Overload these functions.  They are pure virtual in the base class (ATerm).
    //
    virtual casacore::String name() {return casacore::String("EVLA Aperture");};

    virtual void makeFullJones(casacore::ImageInterface<casacore::Complex>& pbImage,
			       const VisBuffer2& vb,
			       casacore::Bool doSquint, casacore::Int& bandID, casacore::Double freqVal);

    virtual void applySky(casacore::ImageInterface<casacore::Float>& outputImages,
			  const VisBuffer2& vb, 
			  const casacore::Bool doSquint=true,
			  const casacore::Int& cfKey=0,
			  const casacore::Int& muellerTerm=0,
			  const casacore::Double freqVal=-1.0);
    virtual void applySky(casacore::ImageInterface<casacore::Complex>& outputImages,
			  const VisBuffer2& vb, 
			  const casacore::Bool doSquint=true,
			  const casacore::Int& cfKey=0,
			  const casacore::Int& muellerTerm=0,
			  const casacore::Double freqVal=-1.0);
    virtual void applySky(casacore::ImageInterface<casacore::Complex>& outImages,
			  const casacore::Double& pa,
			  const casacore::Bool doSquint,
			  const casacore::Int& cfKey,
			  const casacore::Int& muellerTerm,
			  const casacore::Double freqVal=-1.0);
    //Average RR and LL beam 
    //This is the best to avoid considering the squint of the VLA R and L receptors 
    void applyAvgSkyJones(casacore::ImageInterface<casacore::Complex> &outImages);
    //This will generate  the diagonal BeamSkyJones (RR, RL, LR, LL) including the squint at given pa
    // This returns the diagonal Mueller terms PB in the 4 planes  (i.e Mueller term 0, 5, 10 , 15)
    // this is not for heterogenous antennas R and L are assumed to come from exactly similar
    //  antennas
    void applyDiagSkyJones(casacore::ImageInterface<casacore::Complex>& outImages,
                                                 const casacore::Double pa);

    // you have to call the following to set between VLA and EVLA (otherwise default is EVLA)
    void cacheVBInfo(const casacore::String& telescopeName, const casacore::Float& diameter);
    void cacheVBInfo(const VisBuffer2& vb);
    casacore::Int getBandID(const casacore::Double& freq, const casacore::String& telescopeName, const casacore::String& bandName);
    //As this is a specialization for VLA or EVLA
    static casacore::String getVLABandName(const casacore::Double& freq, const casacore::String& tel="EVLA");
    casacore::Int getBandID(const casacore::Double& freq, const casacore::String& bandname="");
    virtual casacore::Vector<casacore::Int> vbRow2CFKeyMap(const VisBuffer2& vb, casacore::Int& nUnique)
    {casacore::Vector<casacore::Int> tmp; tmp.resize(vb.nRows()); tmp=0; nUnique=1; return tmp;}

    virtual void getPolMap(casacore::Vector<casacore::Int>& polMap) {polMap.resize(0);polMap=polMap_p;};

    // For this class, these will be served from the base classs (ATerm.h)
    // virtual casacore::Int getConvSize() {return CONVSIZE;};
    // virtual casacore::Int getOversampling() {return OVERSAMPLING;}
    // virtual casacore::Float getConvWeightSizeFactor() {return CONVWTSIZEFACTOR;};
    // virtual casacore::Float getSupportThreshold() {return THRESHOLD;};

  protected:
    int getVisParams(const VisBuffer2& vb,const casacore::CoordinateSystem& skyCoord=casacore::CoordinateSystem());
    casacore::Bool findSupport(casacore::Array<casacore::Complex>& func, casacore::Float& threshold,casacore::Int& origin, casacore::Int& R);
    casacore::Int getVLABandID(casacore::Double& freq,casacore::String&telescopeName, const casacore::CoordinateSystem& skyCoord=casacore::CoordinateSystem());
    casacore::Int makePBPolnCoords(const VisBuffer2&vb,
			 const casacore::Int& convSize,
			 const casacore::Int& convSampling,
			 const casacore::CoordinateSystem& skyCoord,
			 const casacore::Int& skyNx, const casacore::Int& skyNy,
			 casacore::CoordinateSystem& feedCoord);
    void setApertureParams(ApertureCalcParams& ap,
			   const casacore::Float& Freq, const casacore::Float& pa, 
			   const casacore::Int& bandID,
			   const casacore::IPosition& skyShape,
			   const casacore::Vector<casacore::Double>& uvIncr);

  private:
    casacore::Vector<casacore::Int> polMap_p;
    casacore::Vector<casacore::Int> feedStokes_p;
  };
};
};
#endif
