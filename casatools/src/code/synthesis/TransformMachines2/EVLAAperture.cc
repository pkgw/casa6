// -*- C++ -*-
//# EVLAAperture.cc: Implementation of the EVLAAperture class
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
//#        Internet email: casa-feedback@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$
//
#include <casacore/coordinates/Coordinates/DirectionCoordinate.h>
#include <casacore/coordinates/Coordinates/LinearCoordinate.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/StokesCoordinate.h>
#include <casacore/images/Images/SubImage.h>
#include <casacore/ms/MeasurementSets/MSColumns.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/VisibilityIterator2.h>

#include <synthesis/TransformMachines/SynthesisError.h>
#include <synthesis/TransformMachines2/EVLAAperture.h>
#include <synthesis/TransformMachines2/Utils.h>
#include <synthesis/TransformMachines2/VLACalcIlluminationConvFunc.h>
#include <synthesis/TransformMachines2/WTerm.h>
#include <synthesis/Utilities/FFT2D.h>
//
//---------------------------------------------------------------------
//

using namespace casacore;
namespace casa {
using namespace vi;
using namespace refim;
using namespace SynthesisUtils;

EVLAAperture::EVLAAperture() : AzElAperture(), polMap_p(), feedStokes_p() {
  telescopeName_p = "EVLA";
  Diameter_p = 25.0;
}

EVLAAperture &EVLAAperture::operator=(const EVLAAperture &other) {
  if (this != &other) {
    //	ConvolutionFunction::operator=(other);
    logIO_p = other.logIO_p;
    //	setParams(other.polMap_p_base, other.feedStokes_p);
    setPolMap(other.polMap_p_base);
    telescopeName_p = other.telescopeName_p;
    Diameter_p = other.Diameter_p;
    Nant_p = other.Nant_p;
    HPBW = other.HPBW;
    sigma = other.sigma;
  }
  return *this;
}
Int EVLAAperture::getVLABandID(Double &vbRefFreq, String &telescopeName,
                               const CoordinateSystem &skyCoord) {
  LogIO log_l(LogOrigin("EVLAAperture", "getVLABandID[R&D]"));

  Double refFreq =
      skyCoord.spectralCoordinate(skyCoord.findCoordinate(Coordinate::SPECTRAL))
          .referenceValue()(0);
  //cerr << "getVLABand (Global VB Ref. min, CF Ref.): " << vbRefFreq << " , "
  //     << refFreq << endl;
  if (telescopeName == "VLA") {
    if ((refFreq >= 1.34E9) && (refFreq <= 1.73E9))
      return BeamCalc_VLA_L;
    else if ((refFreq >= 4.5E9) && (refFreq <= 5.0E9))
      return BeamCalc_VLA_C;
    else if ((refFreq >= 8.0E9) && (refFreq <= 8.8E9))
      return BeamCalc_VLA_X;
    else if ((refFreq >= 14.4E9) && (refFreq <= 15.4E9))
      return BeamCalc_VLA_U;
    else if ((refFreq >= 22.0E9) && (refFreq <= 24.0E9))
      return BeamCalc_VLA_K;
    else if ((refFreq >= 40.0E9) && (refFreq <= 50.0E9))
      return BeamCalc_VLA_Q;
    else if ((refFreq >= 100E6) && (refFreq <= 300E6))
      return BeamCalc_VLA_4;
  } else if (telescopeName == "EVLA") {
    if ((refFreq >= 0.9E9) && (refFreq <= 2.1E9))
      return BeamCalc_EVLA_L;
    else if ((refFreq >= 2.0E9) && (refFreq <= 4.0E9))
      return BeamCalc_EVLA_S;
    else if ((refFreq >= 4.0E9) && (refFreq <= 8.0E9))
      return BeamCalc_EVLA_C;
    else if ((refFreq >= 8.0E9) && (refFreq <= 12.0E9))
      return BeamCalc_EVLA_X;
    else if ((refFreq >= 12.0E9) && (refFreq <= 18.0E9))
      return BeamCalc_EVLA_U;
    else if ((refFreq >= 18.0E9) && (refFreq <= 26.5E9))
      return BeamCalc_EVLA_K;
    else if ((refFreq >= 26.5E9) && (refFreq <= 40.8E9))
      return BeamCalc_EVLA_A;
    else if ((refFreq >= 40.0E9) && (refFreq <= 50.0E9))
      return BeamCalc_EVLA_Q;
  }
  ostringstream mesg;
  log_l << telescopeName << "/" << refFreq << "(Hz) combination not recognized."
        << LogIO::EXCEPTION;
  return -1;
}

void EVLAAperture::setApertureParams(ApertureCalcParams &ap, const Float &Freq,
                                     const Float &pa, const Int &bandID,
                                     const IPosition &skyShape,
                                     const Vector<Double> &uvIncr) {
  Double Lambda = C::c / Freq;

  ap.oversamp = 3; 
  ap.pa = pa;
  ap.band = bandID;
  ap.freq = Freq / 1E9;
  ap.nx = skyShape(0);
  ap.ny = skyShape(1);
  ap.dx = abs(uvIncr(0) * Lambda);
  ap.dy = abs(uvIncr(1) * Lambda);
  ap.x0 = -(ap.nx / 2) * ap.dx;
  ap.y0 = -(ap.ny / 2) * ap.dy;
  //cerr << "pa= " << ap.pa << " band " << ap.band << " freq " << ap.freq
  //     << " nx ny " << ap.nx << "  " << ap.ny << " dx dy " << ap.dx << "  "
  //     << ap.dy << endl;
  //
  // If cross-hand pols. are requested, we need to compute both
  // the parallel-hand aperture illuminations.
  //
  // if ((inStokes == Stokes::RL) || (inStokes == Stokes::LR))
  {
    // IPosition apShape(ap.aperture->shape());
    // cerr << "APshape " << apShape << endl;
    // apShape(3)=4;
    //cerr << "APSHPE=" << skyShape << endl;
    
    //ap.aperture->resize(skyShape);
  }
}
void EVLAAperture::cacheVBInfo(const String &telescopeName,
                               const Float &diameter) {
  telescopeName_p = telescopeName;
  Diameter_p = diameter;
}

void EVLAAperture::cacheVBInfo(const VisBuffer2 &vb) {
  const Vector<String> telescopeNames =
      vb.subtableColumns().observation().telescopeName().getColumn();
  for (uInt nt = 0; nt < telescopeNames.nelements(); nt++) {
    if ((telescopeNames(nt) != "VLA") && (telescopeNames(nt) != "EVLA")) {
      String mesg = "We can handle only (E)VLA antennas for now.\n";
      mesg += "Erroneous telescope name = " + telescopeNames(nt) + ".";
      SynthesisError err(mesg);
      throw(err);
    }
    if (telescopeNames(nt) != telescopeNames(0)) {
      String mesg =
          "We do not (yet) handle multiple telescopes for A-Projection!\n";
      mesg += "Not yet a \"priority\"!!";
      SynthesisError err(mesg);
      throw(err);
    }
  }
  telescopeName_p = telescopeNames[0];

  //    MSSpWindowColumns mssp(vb.msColumns().spectralWindow());
  // Freq = vb.msColumns().spectralWindow().refFrequency()(0);
  Diameter_p = 0;
  Nant_p = vb.subtableColumns().antenna().nrow();
  for (Int i = 0; i < Nant_p; i++)
    if (!vb.subtableColumns().antenna().flagRow()(i)) {
      Diameter_p = vb.subtableColumns().antenna().dishDiameter().getColumn()(i);
      break;
    }
  if (Diameter_p == 0) {
    logIO() << LogOrigin("EVLAAperture", "cacheVBInfo")
            << "No valid or finite sized antenna found in the antenna table. "
            << "Assuming diameter = 25m." << LogIO::WARN << LogIO::POST;
    Diameter_p = 25.0;
  }
  cacheVBInfo(telescopeNames[0], Diameter_p);
}

Int EVLAAperture::getBandID(const Double &freq,
                            const String & /*telescopeName*/,
                            const String &bandName) {
  Int bandID = 0;
  if (!isNoOp()) {
    // First #-separated token in bandName_p is the name of the band used
    Vector<String> tokens = SynthesisUtils::parseBandName(bandName);
    Double elfreq = freq;
    if(telescopeName_p=="VLA" && freq < 1.34e9)   //Some unit test data for VLA are below 1.0 GHz
      elfreq = 1.4e9;
    bandID = BeamCalc::Instance()->getBandID(elfreq, telescopeName_p, tokens(0));
  }

  return bandID;
};
String EVLAAperture::getVLABandName(const Double& freq,  const String& telescopeName) {
  double tol = FLT_EPSILON;
  String bandName = "EVLA_L";
  if (telescopeName == "VLA") {
//    if ((freq >= 1.34E9) && (freq <= 1.73E9))   some unit test data goes from 1 to 2 GHz for VLA !

    if((freq >= (9E8-tol)) && (freq <= (2.0E9+tol)))

      bandName = "VLA_L";
    else if ((freq >= (4.5E9-tol)) && (freq <= (5.0E9+tol)))
      bandName = "VLA_C";
    else if ((freq >= (8.0E9-tol)) && (freq <= (8.8E9+tol)))
      bandName = "VLA_X";
    else if ((freq >= (14.4E9-tol)) && (freq <= (15.4E9+tol)))
      bandName = "VLA_U";
    else if ((freq >= (22.0E9-tol)) && (freq <= (24.0E9+tol)))
      bandName = "VLA_K";
    else if ((freq >= (40.0E9-tol)) && (freq <= (50.0E9+tol)))
      bandName = "VLA_Q";
    else if ((freq >= (30E6-tol)) && (freq <= (100E6+tol)))
      bandName = "VLA_4";
    else
      throw(
          AipsError("Don't know VLA band for frequency=" + String::toString(freq)));
  } else if (telescopeName == "EVLA") {
    if (freq > (9e8-tol) && freq <= 2.0e9)
      bandName = "EVLA_L";
    else if (freq > 2.0e9 && freq <= 4.0e9)
      bandName = "EVLA_S";
    else if (freq > 4.0e9 && freq <= 8.0e9)
      bandName = "EVLA_C";
    else if (freq > 8.e9 && freq <= 12.0e9)
      bandName = "EVLA_X";
    else if (freq > 12.0e9 && freq <= 18.0e9)
      bandName = "EVLA_U";
    else if (freq > 18.0e9 && freq <= 26.0e9)
      bandName = "EVLA_K";
    else if (freq > 26.e9 && freq <= 40.0e9)
      bandName = "EVLA_A";
    else if (freq > 40.0e9 && freq <= (50.0e9+tol))
      bandName = "EVLA_Q";
    else
      throw(
          AipsError("Don't know EVLA band for frequency=" + String::toString(freq)));
  } else {
    throw(AipsError("Don't know telescope " + telescopeName));
  }
  return bandName;
}
Int EVLAAperture::getBandID(const Double &freq, const String& bandName) {
  Int bandID = 0;
  if(bandName=="")
    bandName_p = getVLABandName(freq,  telescopeName_p);
  else
    bandName_p=bandName;

  if (!isNoOp()) {
    // First #-separated token in bandName_p is the name of the band used
    // Vector<String> tokens = SynthesisUtils::parseBandName(bandName_p);
    // cerr << "TOKENS " << tokens << endl;
    Double elfreq = freq;
    if (telescopeName_p == "VLA" &&
        freq < 1.34e9) // Some unit test data for VLA are below 1.0 GHz
      elfreq = 1.4e9;
    bandID = BeamCalc::Instance()->getBandID(elfreq, telescopeName_p, bandName_p);
  }

  return bandID;
};
int EVLAAperture::getVisParams(const VisBuffer2 &vb,
                               const CoordinateSystem & /*im*/) {
  throw(AipsError("EVLAAperture::getVisParams() called"));
  Double Freq;
  cacheVBInfo(vb);

  Freq = vb.getFrequency(0, 0);
  Double Lambda = C::c / Freq;
  HPBW = Lambda / (Diameter_p * sqrt(log(2.0)));
  sigma = 1.0 / (HPBW * HPBW);

  return getBandID(Freq);
  // Int bandID=0;
  // if (!isNoOp())
  //   bandID = BeamCalc::Instance()->getBandID(Freq,telescopeName_p);

  // return bandID;
}

Int EVLAAperture::makePBPolnCoords(const VisBuffer2 &vb, const Int &convSize,
                                   const Int &convSampling,
                                   const CoordinateSystem &skyCoord,
                                   const Int &skyNx, const Int & /*skyNy*/,
                                   CoordinateSystem &feedCoord)
//				     Vector<Int>& cfStokes)
{
  feedCoord = skyCoord;
  //
  // Make a two dimensional image to calculate auto-correlation of
  // the ideal illumination pattern. We want this on a fine grid in
  // the UV plane
  //
  Int directionIndex = skyCoord.findCoordinate(Coordinate::DIRECTION);
  //    cout<<"Direction index is "<< directionIndex<<endl;
  AlwaysAssert(directionIndex >= 0, AipsError);
  DirectionCoordinate dc = skyCoord.directionCoordinate(directionIndex);
  Vector<Double> sampling;
  sampling = dc.increment();
  //    cout<<"Image sampling set to : "<<sampling<<endl;
  sampling *= Double(convSampling);
  sampling *= Double(skyNx) / Double(convSize);
  dc.setIncrement(sampling);
  //    cout<<"Resized sampling is : "<<sampling<<endl;

  Vector<Double> unitVec(2);
  unitVec = convSize / 2;
  dc.setReferencePixel(unitVec);

  // Set the reference value to that of the image
  feedCoord.replaceCoordinate(dc, directionIndex);

  //
  // Make an image with circular polarization axis.
  //
  Int NPol = 0, M, N = 0;
  M = polMap_p_base.nelements();
  for (Int i = 0; i < M; i++)
    if (polMap_p_base(i) > -1)
      NPol++;
  Vector<Int> poln(NPol);

  Int index;
  Vector<Int> inStokes;
  index = feedCoord.findCoordinate(Coordinate::STOKES);
  inStokes = feedCoord.stokesCoordinate(index).stokes();
  N = 0;
  try {
    //	cerr << "### " << polMap_p_base << " " << vb.corrType() << endl;
    for (Int i = 0; i < M; i++)
      if (polMap_p_base(i) > -1) {
        poln(N) = vb.correlationTypes()(i);
        N++;
      }
    StokesCoordinate polnCoord(poln);
    Int StokesIndex = feedCoord.findCoordinate(Coordinate::STOKES);
    feedCoord.replaceCoordinate(polnCoord, StokesIndex);
    //	cfStokes = poln;
  } catch (AipsError &x) {
    throw(SynthesisFTMachineError(
        "Likely cause: Discrepancy between the poln. "
        "axis of the data and the image specifications."));
  }

  return NPol;
}

Bool EVLAAperture::findSupport(Array<Complex> &func, Float &threshold,
                               Int &origin, Int &R) {
  Double NSteps;
  Int PixInc = 1;
  Vector<Complex> vals;
  IPosition ndx(4, origin, 0, 0, 0);
  Bool found = false;
  IPosition cfShape = func.shape();
  Int convSize = cfShape(0);
  for (R = convSize / 4; R > 1; R--) {
    NSteps = 90 * R / PixInc; // Check every PixInc pixel along a
    // circle of radious R
    vals.resize((Int)(NSteps + 0.5));
    vals = 0;
    for (Int th = 0; th < NSteps; th++) {
      ndx(0) = (int)(origin + R * sin(2.0 * M_PI * th * PixInc / R));
      ndx(1) = (int)(origin + R * cos(2.0 * M_PI * th * PixInc / R));

      if ((ndx(0) < cfShape(0)) && (ndx(1) < cfShape(1)))
        vals(th) = func(ndx);
    }
    if (max(abs(vals)) > threshold) {
      found = true;
      break;
    }
  }
  return found;
}

void EVLAAperture::makeFullJones(ImageInterface<Complex> &pbImage,
                                 const VisBuffer2 &vb, Bool doSquint,
                                 Int &bandID, Double freqVal) {

  if (!isNoOp()) {
    VLACalcIlluminationConvFunc vlaPB;
    Long cachesize = (HostInfo::memoryTotal(true) / 8) * 1024;
    vlaPB.setMaximumCacheSize(cachesize);
    bandID = getBandID(freqVal);
    vlaPB.makeFullJones(pbImage, vb, doSquint, bandID, freqVal);
  }
}

void EVLAAperture::applyAvgSkyJones(ImageInterface<Complex> &outImage) {
  TempImage<Complex> temp(outImage.shape(), outImage.coordinates());
  temp.setMiscInfo(outImage.miscInfo());
  temp.set(1.0);
  applyDiagSkyJones(temp, 0.0);
  // Taking the RL beam
  IPosition blc(4, 0, 0, 1, 0);
  IPosition trc=outImage.shape()-1;
  trc(2)=1;
  Slicer sl(blc, trc, Slicer::endIsLast);
  SubImage<Complex> rl(temp, sl, false);
  for (Int k = 0; k <4; ++k ) {
   trc[2] = k;
   blc[2] = k;
   sl = Slicer(blc, trc, Slicer::endIsLast);
   SubImage<Complex> outsub(outImage,  sl,  true);
   outsub.copyData(LatticeExpr<Complex>(real(rl)));
  }
  
  
}

void EVLAAperture::applyDiagSkyJones(ImageInterface<Complex> &outImage,
                                     const Double pa) {

  ApertureCalcParams ap;
  Long memtot=HostInfo::memoryFree();
  Double memtobeused= Double(memtot)*1024.0;
  String bandname="";
  //cerr << "diagMisc " << outImage.miscInfo() << endl;
  if(outImage.miscInfo().isDefined("bandname"))
    outImage.miscInfo().get("bandname", bandname);
  //cerr << "diagSky BANDNAME " << bandname << endl;

  ap.apertureptr = make_shared< TempImage<Complex> >(outImage.shape(),  outImage.coordinates(), memtobeused/10.0);  
  ap.aperture = ap.apertureptr.get();
  CoordinateSystem csys = outImage.coordinates();

  Int index = csys.findCoordinate(Coordinate::SPECTRAL);
  Double freqVal = 0.0;
  csys.spectralCoordinate(index).toWorld(freqVal, 0.0);
  //cerr <<  "####FREQVAL " <<  freqVal <<  endl;
  Int bandID = getBandID(freqVal, bandname);
  //cerr << "bandid " << bandID << endl;
  CoordinateSystem uvCoords = refim::VLACalcIlluminationConvFunc::makeUVCoords(
      csys, outImage.shape(), freqVal);
  index = uvCoords.findCoordinate(Coordinate::LINEAR);
  if (index < 0)
    throw(AipsError("Could not get uv-coordinates in applyDiagSkyJones"));
  LinearCoordinate lc = uvCoords.linearCoordinate(index);
  Vector<Double> uvIncr = lc.increment();
  setApertureParams(ap, freqVal, pa, bandID, outImage.shape(), uvIncr);
  //ap.aperture->setCoordinateInfo(uvCoords);
  // Can replace this call with PBMathInterface::applyVP for other type of known beams
  BeamCalc::Instance()->calculateAperture(&ap);
  //cerr <<  "MAx MIN " <<  max(ap.aperture->get()) <<  " " <<  min(ap.aperture->get()) << endl;
  // need to ft to voltagepattern
  FFT2D ft;
  ft.c2cFFT(*(ap.apertureptr));
  // diagonal Mueller beam terms
  // First lets normalize voltage beam
  // In VLACalcIlluminationConvFunc::skyMuller function there is a double conj
  //  first when normalizing tmp ...then when calculation M0 etc
  // So here we are just doing the conj once at multiplication
  // Copying the normalization from skyMuller ...hey it was tested in python must be right
  IPosition mid = outImage.shape()-1;
  mid[0] = outImage.shape()[0]/2;
  mid[1] = outImage.shape()[1]/2;
  Float normalizesq = 0.0;
  for (Int k = 0; k < 4; ++k ) {
    mid[2] = k;
    Complex valmid = ap.aperture->getAt(mid);
    //cerr << " mid " <<  mid <<  " val " <<  valmid <<  endl;
    normalizesq +=  abs(valmid*valmid)/2.0;
  }
  if (normalizesq == 0.0) {
   throw(AipsError("Voltage patterns are 0 at the center")); 
  }
  (ap.aperture)-> copyData(LatticeExpr<Complex>((*(ap.aperture))/sqrt(normalizesq)));
  // now lets calculate the diag Mueller terms remember we are conj the multipliers to match the 
  // double conj that VLACalcIlluminationConvFunc does !
  // For heterogenous arrays that is where we'll multiply the voltage pattern of ant1 with that 
  // of ant2
  
  {
    //  We'll overwrite the 2 middle planes with Jp*conj(Jq) and Jq*conj(Jp) (see skyMuller)
    // because we did not do a conj while normalizing it will be conj(Jp)*Jq and conj(Jq)*Jp 
    //  respectively
    IPosition blc(4, 0, 0, 0, 0);
    IPosition trc=outImage.shape()-1;
    trc(2)=0;
    Slicer sl(blc, trc, Slicer::endIsLast);
    SubImage<Complex> jp(*(ap.aperture), sl, true);
    blc(2) = 3;
    trc(2) = 3;
    sl = Slicer(blc, trc, Slicer::endIsLast);
    SubImage<Complex> jq(*(ap.aperture), sl, true);
    blc(2) = 1;
    trc(2) = 1;
    sl = Slicer(blc, trc, Slicer::endIsLast);
    SubImage<Complex> jpq(*(ap.aperture), sl, true);
    jpq.copyData(LatticeExpr<Complex>(jq*conj(jp)));
    blc(2) = 2;
    trc(2) = 2;
    sl = Slicer(blc, trc, Slicer::endIsLast);
    SubImage<Complex> jqp(*(ap.aperture), sl, true);
    jqp.copyData(LatticeExpr<Complex>(conj(jq)*jp));
    // Now the first and last plane
    jp.copyData(LatticeExpr<Complex>(conj(jp)*jp));
    jq.copyData(LatticeExpr<Complex>(conj(jq)*jq));
  }
  //cerr <<  "SHAPES " <<  ap.aperture->shape() <<  "  " <<  outImage.shape() <<  endl;
  LatticeExpr<Complex> le((*ap.aperture) * outImage);
  outImage.copyData(le);
}

void EVLAAperture::applySky(ImageInterface<Complex> &outImages,
                            const VisBuffer2 &vb, const Bool doSquint,
                            const Int &cfKey, const Int &muellerTerm,
                            const Double freqVal) {
  Double freq = freqVal;
  if (freq < 0.0)
    freq = vb.getFrequency(0, 0);
  //     (void)cfKey;
  //     if (!isNoOp())
  //       {
  // 	VLACalcIlluminationConvFunc vlaPB;
  // 	Long cachesize=(HostInfo::memoryTotal(true)/8)*1024;
  // 	vlaPB.setMaximumCacheSize(cachesize);
  // 	bandID = getBandID(freqVal,telescopeName_p);
  // //	cout<<"EVLAAperture : muellerTerm"<<muellerTerm <<"\n";
  // 	//vlaPB.applyPB(outImages, doSquint,bandID,muellerTerm,freqVal);
  // 	Double pa=getPA(vb);
  // 	vlaPB.applyPB(outImages, pa, doSquint,bandID,muellerTerm,freqVal);

  //       }
  Double pa = getPA(vb);
  //cerr << "PA=" << pa << endl;
  applySky(outImages, pa, doSquint, cfKey, muellerTerm, freq);
}
//
// This one without VB.  Should become the default (and the only one!).
//
void EVLAAperture::applySky(ImageInterface<Complex> &outImages,
                            const Double &pa, const Bool doSquint,
                            const Int &cfKey, const Int &muellerTerm,
                            const Double freqVal) {
  (void)cfKey;
  if (!isNoOp()) {
    VLACalcIlluminationConvFunc vlaPB;
    Long cachesize = (HostInfo::memoryTotal(true) / 8) * 1024;
    vlaPB.setMaximumCacheSize(cachesize);
    Int bandID;
    bandID = getBandID(freqVal);
    // vlaPB.applyPB(outImages, doSquint,bandID,muellerTerm,freqVal);
    Double pa_l = pa; // Due to goofup in making sure complier type checking
                      // does not come in the way!
    //cerr << "PA_L=" << pa_l << " bandID= " << bandID << endl;
    vlaPB.applyPB(outImages, pa_l, doSquint, bandID, muellerTerm, freqVal);
  }
}

void EVLAAperture::applySky(ImageInterface<Float> &outImages,
                            const VisBuffer2 &vb, const Bool doSquint,
                            const Int &cfKey, const Int &muellerTerm,
                            const Double freqVal) {
  (void)cfKey;
  (void)muellerTerm;
  Double freq = freqVal;
  if (freq < 0.0)
    freq = vb.getFrequency(0, 0);
  if (!isNoOp()) {
    VLACalcIlluminationConvFunc vlaPB;
    Long cachesize = (HostInfo::memoryTotal(true) / 8) * 1024;
    vlaPB.setMaximumCacheSize(cachesize);
    Int bandID;
    bandID = getBandID(freq);

    Double pa = getPA(vb);
    vlaPB.applyPB(outImages, pa, bandID, doSquint, freq);
  }
}

}; // namespace casa
