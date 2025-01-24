//# PhaseShiftingTVI.h: This file contains the interface definition of the MSTransformManager class.
//#
//#  CASA - Common Astronomy Software Applications (http://casa.nrao.edu/)
//#  Copyright (C) Associated Universities, Inc. Washington DC, USA 2011, All rights reserved.
//#  Copyright (C) European Southern Observatory, 2011-2024, All rights reserved.
//#
//#  This library is free software; you can redistribute it and/or
//#  modify it under the terms of the GNU Lesser General Public
//#  License as published by the Free software Foundation; either
//#  version 2.1 of the License, or (at your option) any later version.
//#
//#  This library is distributed in the hope that it will be useful,
//#  but WITHOUT ANY WARRANTY, without even the implied warranty of
//#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//#  Lesser General Public License for more details.
//#
//#  You should have received a copy of the GNU Lesser General Public
//#  License along with this library; if not, write to the Free Software
//#  Foundation, Inc., 59 Temple Place, Suite 330, Boston,
//#  MA 02111-1307  USA
//# $Id: $

#ifndef PhaseShiftingTVI_H_
#define PhaseShiftingTVI_H_

#include <unordered_map>

// Base class
#include <mstransform/TVI/FreqAxisTVI.h>

// To handle variant parameters
#include <stdcasa/StdCasa/CasacSupport.h>

// CAS-12706 To access directly the MS columns to get initial
// time and position for wide-field phase shifting algorithm
#include <casacore/ms/MeasurementSets/MSColumns.h>

// CAS-12706 UVWMachine for wide-field phase shifting algorithm
#include <casacore/measures/Measures/UVWMachine.h>
#include <casacore/measures/Measures/MeasFrame.h>


namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE VI - BEGIN

//////////////////////////////////////////////////////////////////////////
// PhaseShiftingTVI class
//////////////////////////////////////////////////////////////////////////

class PhaseShiftingTVI : public FreqAxisTVI
{

public:

	PhaseShiftingTVI(	ViImplementation2 * inputVii,
						const casacore::Record &configuration);

	// Report the the ViImplementation type
	virtual casacore::String ViiType() const { return casacore::String("PhaseShifting( ")+getVii()->ViiType()+" )"; };

	// Navigation methods
	virtual void origin ();
	virtual void next ();

    virtual void visibilityObserved (casacore::Cube<casacore::Complex> & vis) const;
    virtual void visibilityCorrected (casacore::Cube<casacore::Complex> & vis) const;
    virtual void visibilityModel (casacore::Cube<casacore::Complex> & vis) const;
    virtual void uvw (casacore::Matrix<double> & uvw) const;

protected:

    void parseConfiguration(const casacore::Record &configuration);
    void parsePhasecenter(const casacore::Record &configuration);
    void parsePhasecenterDict(const casacore::Record &configuration);
    casacore::MDirection checkPhaseCenterStr(const casacore::String &phasecenter,
                                             const std::string &fieldInfo = "");
    casacore::rownr_t getMaxMSFieldID() const;
    std::pair<bool, casacore::MDirection> findConvertedPhaseCenter() const;
    void initialize();
    void shiftUVWPhases();

	// only used when not using wide-field algorithm (not wideFieldMode_p)
	casacore::Double dx_p, dy_p;

	std::unordered_map<int, casacore::MDirection> phaseCenterSpec_p;
	// CAS-12706 Members wide-field phase shifting algorithm
	bool wideFieldMode_p;

	// Set once from Vii/MS at init time, if wideFieldMode_p
	casacore::MPosition observatoryPosition_p;
	casacore::MEpoch referenceTime_p;
	casacore::String referenceTimeUnits_p;

	// buffer to pass shift from shiftUVWPhases() => TVI uvw/data cols
	casacore::Matrix<casacore::Double> newUVW_p;
	casacore::Vector<casacore::Double> phaseShift_p;
};

//////////////////////////////////////////////////////////////////////////
// PhaseShiftingTVIFactory class
//////////////////////////////////////////////////////////////////////////

class PhaseShiftingTVIFactory : public ViFactory
{

public:

	PhaseShiftingTVIFactory(casacore::Record &configuration,ViImplementation2 *inputVII);

protected:

	vi::ViImplementation2 * createVi (VisibilityIterator2 *) const;
	vi::ViImplementation2 * createVi () const;

	casacore::Record configuration_p;
	ViImplementation2 *inputVii_p;;
};

//////////////////////////////////////////////////////////////////////////
// PhaseShiftingTVILayerFactory class  (for _recursive_ layering)
//////////////////////////////////////////////////////////////////////////

class PhaseShiftingTVILayerFactory : public ViiLayerFactory
{

public:

  PhaseShiftingTVILayerFactory(casacore::Record &configuration);

  virtual ~PhaseShiftingTVILayerFactory() {};

protected:

  
  virtual ViImplementation2 * createInstance(ViImplementation2* vii0) const;

  const casacore::Record configuration_p;

};

//////////////////////////////////////////////////////////////////////////
// PhaseShiftingTransformEngine class
//////////////////////////////////////////////////////////////////////////

template<class T> class PhaseShiftingTransformEngine : public FreqAxisTransformEngine2<T>
{
	using FreqAxisTransformEngine2<T>::inputData_p;
	using FreqAxisTransformEngine2<T>::outputData_p;
	using FreqAxisTransformEngine2<T>::rowIndex_p;
	using FreqAxisTransformEngine2<T>::corrIndex_p;

public:

	PhaseShiftingTransformEngine(casacore::Double dx, casacore::Double dy,
								casacore::Matrix<casacore::Double> *uvw,
								casacore::Vector<casacore::Double> *frequencies,
								DataCubeMap *inputData,
								DataCubeMap *outputData);

	void transform();

	void transformCore(DataCubeMap *inputData,DataCubeMap *outputData);

protected:

	casacore::Double dx_p, dy_p;
	casacore::Matrix<casacore::Double> *uvw_p;
	casacore::Vector<casacore::Double> *frequencies_p;
};

//////////////////////////////////////////////////////////////////////////
// WideFieldPhaseShiftingTransformEngine class
//////////////////////////////////////////////////////////////////////////

template<class T> class WideFieldPhaseShiftingTransformEngine : public FreqAxisTransformEngine2<T>
{
	using FreqAxisTransformEngine2<T>::inputData_p;
	using FreqAxisTransformEngine2<T>::outputData_p;
	using FreqAxisTransformEngine2<T>::rowIndex_p;
	using FreqAxisTransformEngine2<T>::corrIndex_p;

public:

	WideFieldPhaseShiftingTransformEngine(const casacore::Vector<casacore::Double> &phaseShift,
								casacore::Matrix<casacore::Double> *uvw,
								casacore::Vector<casacore::Double> *frequencies,
								DataCubeMap *inputData,
								DataCubeMap *outputData);

	void transform();

	void transformCore(DataCubeMap *inputData,DataCubeMap *outputData);

protected:

	casacore::Vector<casacore::Double> phaseShift_p;
	casacore::Matrix<casacore::Double> *uvw_p;
	casacore::Vector<casacore::Double> *frequencies_p;
};



} //# NAMESPACE VI - END

} //# NAMESPACE CASA - END

#endif /* PhaseShiftingTVI_H_ */

