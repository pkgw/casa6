//# PhaseShiftingTVI.h: This file contains the implementation of the PhaseShiftingTVI class.
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

#include <mstransform/TVI/PhaseShiftingTVI.h>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN

namespace vi { //# NAMESPACE VI - BEGIN

//////////////////////////////////////////////////////////////////////////
// PhaseShiftingTVI class
//////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
PhaseShiftingTVI::PhaseShiftingTVI(	ViImplementation2 * inputVii,
								const Record &configuration):
								FreqAxisTVI (inputVii)
{
	dx_p = 0;
	dy_p = 0;

	// CAS-12706 Zero initialization for wide-field phase shifting algorithm
	wideFieldMode_p = false;

	// Parse and check configuration parameters
	// Note: if a constructor finishes by throwing an exception, the memory
	// associated with the object itself is cleaned up — there is no memory leak.
	parseConfiguration(configuration);

	initialize();

	return;
}

/**
 * Get max valid FIELD IDs for this iterator. MSv2 uses the FIELD
 * table row index as FIELD ID.
 */
rownr_t PhaseShiftingTVI::getMaxMSFieldID() const
{
  const auto &fieldsTable = getVii()->fieldSubtablecols();
  return fieldsTable.nrow() - 1;
}

/**
 * Parses the phase center parameter, whether given as a single string or
 * as a dictionary/record of per-field strings.
 * Populates the phaseCenterSpec_p map.
 *
 * @param config TVI configuration object
 */
void PhaseShiftingTVI::parsePhasecenter(const Record &config)
{
  auto exists = config.fieldNumber ("phasecenter");
  if (exists < 0) {
    return;
  }

  // phasecenter can be given as a string or as a dict (per-field centers)
  auto isStr = false;
  casacore::String phaseCenterStr;
  try {
    config.get(exists, phaseCenterStr);
    isStr = true;
    logger_p << LogIO::NORMAL << LogOrigin("PhaseShiftingTVI", __FUNCTION__)
	     << "Interpreting phasecenter as a single string, to be applied"
	     " to all input (selected) fields."<< LogIO::POST;
  } catch(const AipsError &exc) {
    // ignore "RecordRep::get_pointer - incorrect data type String used...
    //   for field phasecenter with type Record" and similar, try record type
    logger_p << LogIO::NORMAL << LogOrigin("PhaseShiftingTVI", __FUNCTION__)
	     << "Interpreting phasecenter as a dictionary of field->center."
	     << "Using following phase centers:" << LogIO::POST;
  }

  if (isStr) {
    const auto &phaseCenterDir = checkPhaseCenterStr(phaseCenterStr);
    phaseCenterSpec_p.insert({-1, phaseCenterDir});
  } else {
    parsePhasecenterDict(config);
  }
}

/**
 * Parses the phase center strings in the items of a dict/record.
 * Populates the phaseCenterSpec_p map.
 *
 * @param config TVI configuration object
 */
void PhaseShiftingTVI::parsePhasecenterDict(const Record &config)
{
  const auto &inPhasecenter = config.asRecord("phasecenter");
  if (inPhasecenter.empty()) {
    throw AipsError("The dictionary 'phasecenter' is empty.");
  }

  const auto maxMSField = getMaxMSFieldID();
  std::set<unsigned int> fieldsSeen;
  // Go through items in the input dict/record
  for (unsigned int rid=0; rid < inPhasecenter.nfields(); ++rid) {
    const std::string fieldStr = inPhasecenter.name(RecordFieldId(rid));
    const auto fid = std::stoi(fieldStr);
    if (fid < 0 || static_cast<unsigned int>(fid) > maxMSField) {
      throw AipsError("Wrong field ID given: " + std::to_string(fid) +
		      ". This MeasurementSet has field IDs between 0 and " +
		      std::to_string(maxMSField));
    }
    if (not fieldsSeen.insert(fid).second) {
      throw AipsError("Field " + std::to_string(fid) + " is given multiple times");
    }

    std::string center;
    try {
      center = inPhasecenter.asString(RecordFieldId(rid));
    } catch (const AipsError &exc) {
      throw AipsError("For field " + std::to_string(fid) + ", cannot interpret "
		      "phasecenter value as a string: " +
		      std::string(exc.getMesg()));
    }
    const auto strFieldID = std::to_string(fid);
    const auto &phaseCenterDir = checkPhaseCenterStr(center, strFieldID);
    phaseCenterSpec_p.insert({fid, phaseCenterDir});
  }
}

/**
 * Checks that a phase center string can be correctly converted to a
 * MDirection object, and checks that the ref frame is supported.
 *
 * @param phasecenter center string as given in the 'phasecenter' param
 * @param fieldInfo more info to print about the field
 *
 * @returns direction as an MDirection
 */
MDirection PhaseShiftingTVI::checkPhaseCenterStr(const String &phasecenter,
						 const string &fieldInfo)
{
  // casaMDirection requires a variant
  casac::variant phaseCenterVar(phasecenter);
  casacore::MDirection phaseCenterDir;
  if(!casaMDirection(phaseCenterVar, phaseCenterDir)) {
    throw AipsError("Cannot interpret phase center string as a direction object: "
		    + phasecenter);
  } else {
    const auto myFrame = phaseCenterDir.getRefString();
    MDirection::Types mdtype;
    MDirection::getType(mdtype, myFrame);
    ThrowIf(
	    mdtype == MDirection::HADEC || mdtype == MDirection::AZEL
	    || mdtype == MDirection::AZELSW || mdtype == MDirection::AZELNE
	    || mdtype == MDirection::AZELGEO || mdtype == MDirection::AZELSWGEO
	    || mdtype == MDirection::MECLIPTIC || mdtype == MDirection::TECLIPTIC
	    || mdtype == MDirection::TOPO,
	    myFrame + " is a time dependent reference frame and so is not supported"
	    );
    ThrowIf(
	    mdtype == MDirection::MERCURY || mdtype == MDirection::VENUS
	    || mdtype == MDirection::MARS || mdtype == MDirection::JUPITER
	    || mdtype == MDirection::SATURN || mdtype == MDirection::URANUS
	    || mdtype == MDirection::NEPTUNE || mdtype == MDirection::PLUTO
	    || mdtype == MDirection::SUN || mdtype == MDirection::MOON
	    || mdtype == MDirection::COMET,
	    myFrame + " denotes an ephemeris object and so is not supported"
	    );

    auto msg = "Phase center '" + std::string(phasecenter) +
      "' successfully parsed";
    if (not fieldInfo.empty()) {
      msg +=  " for field " + fieldInfo;
    }
    logger_p << LogIO::NORMAL << LogOrigin("PhaseShiftingTVI", __FUNCTION__)
	     << msg << LogIO::POST;
  }

  return phaseCenterDir;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void PhaseShiftingTVI::parseConfiguration(const Record &configuration)
{
  int exists = -1;

  exists = -1;
  exists = configuration.fieldNumber ("XpcOffset");
  if (exists >= 0)
  {
    configuration.get (exists, dx_p);
  }

  exists = -1;
  exists = configuration.fieldNumber ("YpcOffset");
  if (exists >= 0)
  {
    configuration.get (exists, dy_p);
  }

  if (abs(dx_p) > 0 or abs(dy_p) > 0)
  {
    logger_p << LogIO::NORMAL << LogOrigin("PhaseShiftingTVI", __FUNCTION__)
             << "Phase shift is dx="<< dx_p << " dy=" << dy_p << LogIO::POST;
  }

  // CAS-12706 Add support for shifting across large offset/angles
  parsePhasecenter(configuration);
  wideFieldMode_p = true;
}

/**
 * Finds the phase center (converted as MDirection object) for the current
 * field, converting the ref frame if needed to match the current Vis Buffer
 * ref frame.
 *
 * @returns whether any shift should be applied (if an output phasecenter is
 *          defined for this field + phase center to shift to, for current
 *          field
 */
std::pair<bool, MDirection>
PhaseShiftingTVI::findConvertedPhaseCenter() const
{
  const auto *vb = getVii()->getVisBuffer();

  auto centerIt = phaseCenterSpec_p.find(-1);
  if (centerIt == phaseCenterSpec_p.end()) {
    auto fieldID = vb->fieldId()[0];
    centerIt = phaseCenterSpec_p.find(fieldID);
    if (centerIt == phaseCenterSpec_p.end()) {
      return {false, MDirection()};
    }
  }
  auto convertedPhaseCenter = centerIt->second;

  if (convertedPhaseCenter.getRefString() != vb->phaseCenter().getRefString()) {
    MDirection::Types mdtype;
    MDirection::getType(mdtype, vb->phaseCenter().getRefString());
    convertedPhaseCenter = MDirection::Convert(convertedPhaseCenter, mdtype)();
  }

  return {true, convertedPhaseCenter};
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void PhaseShiftingTVI::initialize()
{
  // Populate nchan input-output maps
 Int spw;
 uInt spw_idx = 0;
 map<Int,vector<Int> >::iterator iter;
 for (iter=spwInpChanIdxMap_p.begin();iter!=spwInpChanIdxMap_p.end();iter++)
 {
   spw = iter->first;
   spwOutChanNumMap_p[spw] = spwInpChanIdxMap_p[spw].size();

   spw_idx++;
 }

 // CAS-12706 Add support for shifting across large offset/angles
 // Access observatory position and observation start (reference) time.
 if (wideFieldMode_p)
 {
   const auto selectedInputMsCols = std::make_unique<MSColumns>(getVii()->ms());
   observatoryPosition_p = selectedInputMsCols->antenna().positionMeas()(0);
   referenceTime_p  = selectedInputMsCols->timeMeas()(0);
   referenceTimeUnits_p = selectedInputMsCols->timeQuant()(0).getUnit();
 }
}

// -----------------------------------------------------------------------
// CAS-12706 Recalculate UVW and Phases according to wide-field phase shifting
// algorithm. This method should do the same as ComponentFTMachine::rotateUVW
// -----------------------------------------------------------------------
void PhaseShiftingTVI::shiftUVWPhases()
{
  // Get input VisBuffer
  VisBuffer2 *vb = getVii()->getVisBuffer();

  bool doShift = false;
  MDirection convertedPhaseCenter;
  std::tie(doShift, convertedPhaseCenter) = findConvertedPhaseCenter();
  if (not doShift) {
      phaseShift_p.resize(0, false);
      newUVW_p = vb->uvw();
      return;
  }

  // Initialize epoch corresponding to current buffer
  // with time reference to the first row in the MS
  MEpoch epoch(Quantity(vb->time()(0),referenceTimeUnits_p),referenceTime_p.getRef());
  MeasFrame refFrame(epoch,observatoryPosition_p);
  UVWMachine uvwMachine(convertedPhaseCenter, vb->phaseCenter(), refFrame,false,false);
  // Initialize phase array and uvw matrix
  phaseShift_p.resize(vb->nRows(),false);
  newUVW_p.resize(vb->uvw().shape(),false);

  // Obtain phase shift and new uvw coordinates
  Vector<Double> dummy(3,0.0);
  double phase2radPerHz = -2.0 * C::pi / C::c;
  for (rownr_t row=0; row<vb->nRows(); row++)
  {
    // Copy current uvw coordinates so that they are not modified
    // Note: Columns in uvw correspond to rows in the main table/VisBuffer!
    dummy = vb->uvw().column(row);

    // Have to change (u,v,w) to (-u,-v,w) because is the convention used by uvwMachine
    dummy(0) = -1*dummy(0);
    dummy(1) = -1*dummy(1);

    // Transform uvw coordinates and obtain corresponding phase shift
    uvwMachine.convertUVW(phaseShift_p(row), dummy);

    // Store new uvw coordinates
    // Have to change back (-u,-v,w) to (u,v,w) because is the convention used by the MS
    dummy(0) = -1*dummy(0);
    dummy(1) = -1*dummy(1);
    newUVW_p.column(row) = dummy;
    // Convert phase shift to radian/Hz
    phaseShift_p(row) = phase2radPerHz*phaseShift_p(row);
  }
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void PhaseShiftingTVI::origin()
{
  // Drive underlying ViImplementation2
 getVii()->origin();

 // CAS-12706 Add support for shifting across large offset/angles
 if (wideFieldMode_p) {
   shiftUVWPhases();
 }

 // Define the shapes in the VB2, patch provided by cgarcia in CAS-12706
 configureShapes();

 // Synchronize own VisBuffer
 configureNewSubchunk();
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void PhaseShiftingTVI::next()
{
  // Drive underlying ViImplementation2
  getVii()->next();

  // CAS-12706 Add support for shifting across large offset/angles
  if (wideFieldMode_p) {
    shiftUVWPhases();
  }

  // Define the shapes in the VB2, patch provided by cgarcia in CAS-12706
  configureShapes();

  // Synchronize own VisBuffer
  configureNewSubchunk();
}


// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void PhaseShiftingTVI::visibilityObserved (Cube<Complex> & vis) const
{
  // Get input VisBuffer
  VisBuffer2 *vb = getVii()->getVisBuffer();
  Matrix<Double> uvw = vb->uvw();
  Vector<Double> frequencies = vb->getFrequencies(0);

  // Reshape output data before passing it to the DataCubeHolder
  vis.resize(getVisBuffer()->getShape(),false);

  // Gather input data
  DataCubeMap inputData;
  DataCubeHolder<Complex> inputVisCubeHolder(vb->visCube());
  inputData.add(MS::DATA,inputVisCubeHolder);

  // Gather output data
  DataCubeMap outputData;
  DataCubeHolder<Complex> outputVisCubeHolder(vis);
  outputData.add(MS::DATA,outputVisCubeHolder);

  if (wideFieldMode_p)
  {
    // Configure Transformation Engine
    WideFieldPhaseShiftingTransformEngine<Complex> transformer(phaseShift_p,&uvw,&frequencies,&inputData,&outputData);

    // Transform data
    transformFreqAxis2(vb->getShape(),transformer);
  }
  else
  {
    // Configure Transformation Engine
    PhaseShiftingTransformEngine<Complex> transformer(dx_p,dy_p,&uvw,&frequencies,&inputData,&outputData);

    // Transform data
    transformFreqAxis2(vb->getShape(),transformer);
  }
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void PhaseShiftingTVI::visibilityCorrected (Cube<Complex> & vis) const
{
  // Get input VisBuffer
 VisBuffer2 *vb = getVii()->getVisBuffer();
 Matrix<Double> uvw = vb->uvw();
 Vector<Double> frequencies = vb->getFrequencies(0);

 // Reshape output data before passing it to the DataCubeHolder
 vis.resize(getVisBuffer()->getShape(),false);

 // Gather input data
 DataCubeMap inputData;
 DataCubeHolder<Complex> inputVisCubeHolder(vb->visCubeCorrected());
 inputData.add(MS::DATA,inputVisCubeHolder);

 // Gather output data
 DataCubeMap outputData;
 DataCubeHolder<Complex> outputVisCubeHolder(vis);
 outputData.add(MS::DATA,outputVisCubeHolder);

 if (wideFieldMode_p)
 {
   // Configure Transformation Engine
   WideFieldPhaseShiftingTransformEngine<Complex> transformer(phaseShift_p,&uvw,&frequencies,&inputData,&outputData);

   // Transform data
   transformFreqAxis2(vb->getShape(),transformer);
 }
 else
 {
   // Configure Transformation Engine
   PhaseShiftingTransformEngine<Complex> transformer(dx_p,dy_p,&uvw,&frequencies,&inputData,&outputData);

   // Transform data
   transformFreqAxis2(vb->getShape(),transformer);
 }
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void PhaseShiftingTVI::visibilityModel (Cube<Complex> & vis) const
{
  // Get input VisBuffer
  VisBuffer2 *vb = getVii()->getVisBuffer();
  Matrix<Double> uvw = vb->uvw();
  Vector<Double> frequencies = vb->getFrequencies(0);

  // Reshape output data before passing it to the DataCubeHolder
  vis.resize(getVisBuffer()->getShape(),false);

  // Gather input data
  DataCubeMap inputData;
  DataCubeHolder<Complex> inputVisCubeHolder(vb->visCubeModel());
  inputData.add(MS::DATA,inputVisCubeHolder);

  // Gather output data
  DataCubeMap outputData;
  DataCubeHolder<Complex> outputVisCubeHolder(vis);
  outputData.add(MS::DATA,outputVisCubeHolder);

  if (wideFieldMode_p)
  {
    // Configure Transformation Engine
    WideFieldPhaseShiftingTransformEngine<Complex> transformer(phaseShift_p,&uvw,&frequencies,&inputData,&outputData);

    // Transform data
    transformFreqAxis2(vb->getShape(),transformer);
  }
  else
  {
      // Configure Transformation Engine
      PhaseShiftingTransformEngine<Complex> transformer(dx_p,dy_p,&uvw,&frequencies,&inputData,&outputData);

      // Transform data
      transformFreqAxis2(vb->getShape(),transformer);
  }
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
void PhaseShiftingTVI::uvw (casacore::Matrix<double> & uvw) const
{
  if (wideFieldMode_p) {
    uvw.resize(newUVW_p.shape(),false);
    uvw = newUVW_p;
  }
  else {
    getVii()->uvw (uvw);
  }
}

//////////////////////////////////////////////////////////////////////////
// PhaseShiftingTVIFactory class
//////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
PhaseShiftingTVIFactory::PhaseShiftingTVIFactory (	Record &configuration,
													ViImplementation2 *inputVii)
{
  inputVii_p = inputVii;
  configuration_p = configuration;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
vi::ViImplementation2 * PhaseShiftingTVIFactory::createVi(VisibilityIterator2 *) const
{
  return new PhaseShiftingTVI(inputVii_p,configuration_p);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
vi::ViImplementation2 * PhaseShiftingTVIFactory::createVi() const
{
  return new PhaseShiftingTVI(inputVii_p,configuration_p);
}

//////////////////////////////////////////////////////////////////////////
// PhaseShiftingTVILayerFactory class
//////////////////////////////////////////////////////////////////////////

PhaseShiftingTVILayerFactory::PhaseShiftingTVILayerFactory(Record &configuration) :
  ViiLayerFactory(),
  configuration_p(configuration)
{}

ViImplementation2* 
PhaseShiftingTVILayerFactory::createInstance(ViImplementation2* vii0) const 
{
  // Make the PhaseShiftingTVi2, using supplied ViImplementation2, and return it
  ViImplementation2 *vii = new PhaseShiftingTVI(vii0,configuration_p);
  return vii;
}

//////////////////////////////////////////////////////////////////////////
// PhaseShiftingTransformEngine class
//////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
template<class T> PhaseShiftingTransformEngine<T>::PhaseShiftingTransformEngine(
															Double dx, Double dy,
															Matrix<Double> *uvw,
															Vector<Double> *frequencies,
															DataCubeMap *inputData,
															DataCubeMap *outputData):
															FreqAxisTransformEngine2<T>(inputData,outputData)
{
  uvw_p = uvw;
  frequencies_p = frequencies;

  // Offsets in radians (input is arcsec)
  dx_p = dx*(C::pi / 180.0 / 3600.0);
  dy_p = dy*(C::pi / 180.0 / 3600.0);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
template<class T> void PhaseShiftingTransformEngine<T>::transform(	)
{
   transformCore(inputData_p,outputData_p);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
template<class T> void PhaseShiftingTransformEngine<T>::transformCore(	DataCubeMap *inputData,
																		DataCubeMap *outputData)
{
  // Get input/output data
  Vector<T> &inputVector = inputData->getVector<T>(MS::DATA);
  Vector<T> &outputVector = outputData->getVector<T>(MS::DATA);

  // Extra path as fraction of U and V in m
  Double phase = dx_p*(*uvw_p)(0,rowIndex_p) + dy_p*(*uvw_p)(1,rowIndex_p);

  // In radian/Hz
  phase *= -2.0 * C::pi / C::c;

  // Main loop
  Double phase_i;
  Complex factor;
  for (uInt chan_i=0;chan_i<inputVector.size();chan_i++)
  {
    phase_i = phase * (*frequencies_p)(chan_i);
    factor = Complex(cos(phase_i), sin(phase_i));
    outputVector(chan_i) = factor*inputVector(chan_i);
  }
}

//////////////////////////////////////////////////////////////////////////
// WideFieldPhaseShiftingTransformEngine class
//////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
template<class T> WideFieldPhaseShiftingTransformEngine<T>::WideFieldPhaseShiftingTransformEngine(
															const Vector<Double> &phaseShift,
															Matrix<Double> *uvw,
															Vector<Double> *frequencies,
															DataCubeMap *inputData,
															DataCubeMap *outputData):
															FreqAxisTransformEngine2<T>(inputData,outputData)
{
  uvw_p = uvw;
  frequencies_p = frequencies;
  phaseShift_p = phaseShift;
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
template<class T> void WideFieldPhaseShiftingTransformEngine<T>::transform(	)
{
  transformCore(inputData_p,outputData_p);
}

// -----------------------------------------------------------------------
//
// -----------------------------------------------------------------------
template<class T> void WideFieldPhaseShiftingTransformEngine<T>::transformCore(	DataCubeMap *inputData,
																		DataCubeMap *outputData)
{
  // Get input/output data
 Vector<T> &inputVector = inputData->getVector<T>(MS::DATA);
 Vector<T> &outputVector = outputData->getVector<T>(MS::DATA);

 if (phaseShift_p.shape() == 0) {
   // no shift, bypass data as 'passthrough' field
   outputVector = inputVector;
   return;
 }

 // Main loop
 Double phase_i;
 Complex factor;
 for (uInt chan_i=0;chan_i<inputVector.size();chan_i++)
 {
   phase_i = phaseShift_p(rowIndex_p) * (*frequencies_p)(chan_i);
   factor = Complex(cos(phase_i), sin(phase_i));
   outputVector(chan_i) = factor*inputVector(chan_i);
 }
}


} //# NAMESPACE VI - END

} //# NAMESPACE CASA - END


