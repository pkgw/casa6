//# AWConvFuncHolser.cc: Implementation for Holding class for AW convolution functions
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
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/VisibilityIterator2.h>
#include <synthesis/TransformMachines2/AWConvFunc.h>
#include <synthesis/TransformMachines2/EVLAAperture.h>
#include <synthesis/TransformMachines2/AWConvFuncHolder.h>

#include <casacore/coordinates/Coordinates/TabularCoordinate.h>
#include <iomanip>
#include <casacore/casa/OS/Timer.h>
namespace casa {  //# CASA namespace
namespace refim { //# namespace refactor imaging

using namespace casacore;
using namespace casa;
using namespace casa::refim;
using namespace std;

AWConvFuncHolder::AWConvFuncHolder(const CoordinateSystem& csys, const int nx, const int ny, const bool dosquint, const double paInc,  const String& obs,  const int& oversamp): painc_p(paInc),  dosquint_p(dosquint), outcsys_p(csys),  nx_p(nx),  ny_p(ny),  oversamp_p(oversamp), isSingleField_p(false) {
  
  convFunc_p.resize();
  wgtConvFunc_p.resize();
  polVals_p.resize();
  freqVals_p.resize();
  wVals_p.resize();
  paVals_p.resize();
  antpairVals_p.resize();
  rowAxisWVals_p.resize();
  rowAxisPAVals_p.resize();
  rowAxisAntennaPair_p.resize();
  convSizes_p.resize();
  convSupport_p.resize();
  convSizesHPG_p.resize();
  convSupportHPG_p.resize();
  aterm_p = std::make_shared<refim::EVLAAperture>();
  aterm_p->cacheVBInfo(obs,  25.0);
  
}
AWConvFuncHolder::AWConvFuncHolder(const AWConvFuncHolder& other) {
 operator = (other);
}
AWConvFuncHolder& AWConvFuncHolder::operator=(const AWConvFuncHolder& other) {
  if (this != &other) {
    nx_p = other.nx_p;
    ny_p = other.ny_p;
    oversamp_p = other.oversamp_p;
    outcsys_p = other.outcsys_p;
    calcNpix_p = other.calcNpix_p;
    calcCsys_p = other.calcCsys_p;
    dosquint_p = other.dosquint_p;
    painc_p = other.painc_p;
    convFunc_p.resize();
    convFunc_p = other.convFunc_p;
    wgtConvFunc_p.resize();
    wgtConvFunc_p = other.wgtConvFunc_p;
    polVals_p.resize();
    polVals_p = other.polVals_p;
    freqVals_p.resize();
    freqVals_p = other.freqVals_p;
    wVals_p.resize();
    wVals_p = other.wVals_p;
    paVals_p.resize();
    paVals_p = other.paVals_p;
    antpairVals_p.resize();
    antpairVals_p = other.antpairVals_p;
    rowAxisWVals_p.resize();
    rowAxisWVals_p = other.rowAxisWVals_p;
    rowAxisPAVals_p.resize();
    rowAxisPAVals_p = other.rowAxisPAVals_p;
    rowAxisAntennaPair_p.resize();
    rowAxisAntennaPair_p = other.rowAxisAntennaPair_p;
    convSizes_p.resize();
    convSizesHPG_p.resize();
    convSizes_p = other.convSizes_p;
    convSizesHPG_p = other.convSizesHPG_p;

    convSupport_p.resize();
    convSupport_p.resize();
    convSupportHPG_p.resize();
    convSupport_p = other.convSupport_p;
    convSupportHPG_p = other.convSupportHPG_p;
    aterm_p = other.aterm_p;
    isSingleField_p = other.isSingleField_p;
  }
  return *this;
}
bool AWConvFuncHolder::addConvFunc(const casacore::Vector<casacore::Double>& freqs, const casacore::Vector<Double>& wVals,  const casacore::Double& paMax) {
  
  Vector<Double> freqsToCalc;
 if (freqVals_p.nelements() == 0) {
   freqVals_p = freqs;
   freqsToCalc = freqs;
   
 }
 else{
  //append new freqs freqVals_p;
  // assign missing freqs to freqsToCalc
  freqsToCalc =freqs;  
 }
 if (wVals_p.nelements() == 0) {
   // make sure first wval is 0;
  wVals_p = wVals;
  
 }
 else if (wVals_p.nelements() !=  wVals.nelements())
   throw(AipsError("Cannot change W terms length right now"));
 if (!dosquint_p) {
   paVals_p.resize(1);
   paVals_p[0] = 0.0;
 }
 else{
   cerr << "paMax " << paMax << " painc " << painc_p << endl;
   Vector<Double> pavals(int(std::ceil(2 * paMax / painc_p)));
   // setting pavals from -paMax to paMax
   for (uint k = 0; k < pavals.nelements(); ++k)
     pavals[k] = double(k) * painc_p - paMax;
   if (paVals_p.nelements() == 0)
     paVals_p = pavals;
   else if ((paVals_p.nelements()) != pavals.nelements())
     throw(AipsError("Cannot change number of PA's in between"));

 }
  std::shared_ptr<refim::WPConvFunc>wptr;
  AWConvFunc a(aterm_p, wptr);
  //Array<Complex> aWConv;
  //Array<Complex> aWwtconv;
  //Matrix<Int> awSupport;
  calcCsys_p = outcsys_p;
  calcNpix_p = min(nx_p,  ny_p);
  //cerr << "PAVALS " <<  paVals_p <<  " dosquint " << dosquint_p <<  endl;
  //cerr << "FREQS " << freqsToCalc << endl;

  for (uint k = 0; k < paVals_p.nelements(); ++k) {
    calcNpix_p = min(nx_p,  ny_p);
    calcCsys_p = outcsys_p;
    Array<Complex> aWConv;
    Array<Complex> aWwtconv;
    Matrix<Int> awSupport;
    a.makeAWConvFunc(aWConv, aWwtconv, calcCsys_p, awSupport, calcNpix_p,
                     freqsToCalc, wVals_p, dosquint_p, paVals_p[k],
                     isSingleField_p);
    //cerr << "######MAX awsupp " << max(awSupport) << endl;
    int startrow = k * wVals_p.nelements();
    appendConvFuncs(aWConv, aWwtconv, awSupport, freqsToCalc, paVals_p[k], startrow);
    // Let's resize for all paVals as resizing is costly
    if (k == 0 && paVals_p.nelements() > 1) {
      IPosition shp = convFunc_p.shape();
      shp[4] = wVals.nelements() * paVals_p.nelements();
      Double memAmt = Double(shp.product()) * 16.0;
      Double memAvl = Double(HostInfo::memoryFree()) * 1024.0;
      //cerr << "Memory needed " << memAmt << " available " << memAvl << endl;
      if (memAmt > 0.8 * memAvl) {
        throw(AipsError(
            "Not enough memory to hold all AW with squint correction convolution functions "));
      }
      convFunc_p.resize(shp, true);
      wgtConvFunc_p.resize(shp, true);
    }

  }

 return true;;
}
void AWConvFuncHolder::appendConvFuncs(const Array<Complex>& awConv,  const Array<Complex>& aWwtConv,  const Matrix<Int>& awsupport,  const Vector<Double>& newFreqs,  const Double paVal, const int startrow) {
  Int nfreqs = awConv.shape()[3];
  // Make sure the polVals are in the stokes used in making convfun
  polVals_p.resize(4);
  // for now antpairVals are for all pair of possible antenna combination
  antpairVals_p.resize(1);
  antpairVals_p[0] = std::pair<int,  int>(-1,  -1);
  
  convertArray(polVals_p,  calcCsys_p.stokesCoordinate(1).stokes());
  IPosition blc = freqVals_p.shape();
  IPosition trc(1, blc[0]+nfreqs-1);
  if (freqVals_p.nelements() !=  0) {
   if (!allEQ(freqVals_p,  newFreqs))
     throw(AipsError("Not implemented appending different freqs yet to conv"));
  }
  else{
    freqVals_p.resize(freqVals_p.nelements()+nfreqs,  True);
    freqVals_p(blc, trc) = newFreqs;
  }
  Int indPA = -1;
  for (uint k = 0; k < paVals_p.nelements(); ++k) {
   if (fabs(paVal-paVals_p[k]) < 1e-4)
     indPA = k;
  }
  if (indPA <0) {
    paVals_p.resize(paVals_p.nelements()+1,  True);
    indPA =  paVals_p.nelements()-1;
   paVals_p[indPA] = paVal;
   
  }
  blc[0] = rowAxisPAVals_p.nelements();
  trc[0] = blc[0]+wVals_p.nelements()-1;
  rowAxisPAVals_p.resize(trc[0]+1,  True);
  rowAxisPAVals_p(blc,  trc).set(indPA);
  rowAxisWVals_p.resize(trc[0]+1,  True);
  Vector<Int> waxis(wVals_p.nelements());
  indgen(waxis);
  rowAxisWVals_p(blc, trc) = waxis;
  //cerr << "rowAxisWVals " << rowAxisWVals_p << endl;
  rowAxisAntennaPair_p.resize(trc[0] + 1, True);
  rowAxisAntennaPair_p(blc, trc).set(0);
  /// Let us rescale convolution function to match nx, ny and incr of image that makes 
  /// uvgrid
  Float factorX=fabs(calcCsys_p.increment()(0)/outcsys_p.increment()(0));
  Float factorY=fabs(calcCsys_p.increment()(1)/outcsys_p.increment()(1));
  //cerr <<  "####Factor " <<  factorX <<  "   " <<  factorY <<  endl;
  factorX = Float(nx_p) *Float(oversamp_p)/Float(calcNpix_p)/factorX;
  factorY = Float(ny_p) *Float(oversamp_p)/Float(calcNpix_p)/factorY;
  
  //cerr <<  "factors " <<  factorX <<  "   " <<  factorY <<  "nx,  ny" <<  nx_p << "   " << ny_p << " calcNpix " << calcNpix_p << " oversamp " << oversamp_p << endl;
  MathUtils m;
  Array<Complex> newAWConv;
  Array<Complex> newWtConv;
  // For small images or factor less than 1.0  use linear interpolation
  //cerr << "Shape before " << awConv.shape() << endl;
  // factor/oversamp is ratio of im fov  to pb fov
  if ((factorX / oversamp_p) < 1.0 || (factorY / oversamp_p) < 1.0 ||
      nx_p < 200 || ny_p < 200) {
    newAWConv = m.resample(awConv, factorX, factorY);
    newWtConv = m.resample(aWwtConv, factorX, factorY);


  } else {
    newAWConv = m.resampleViaFFT(awConv, factorX, factorY);
    newWtConv = m.resampleViaFFT(aWwtConv,  factorX, factorY);
  }
  //cerr << "Shape after " << newAWConv.shape() << endl;
  Float correcfac =
      float(awConv.shape()(0) * awConv.shape()(1) * oversamp_p * oversamp_p) /
      float(newAWConv.shape()(0) * newAWConv.shape()(1));

  newAWConv *= correcfac;
  newWtConv *= correcfac;

  /*{ 
      ////TESTOO
      IPosition elshp = newAWConv.shape().getFirst(4);
      IPosition elblc(5, 0);
      
      IPosition eltrc = newAWConv.shape()-1;
      elblc[4] = eltrc[4];
      PagedImage<Complex> lastplane(elshp,  calcCsys_p,  "MOOBOO");
      lastplane.put(newAWConv(elblc, eltrc).nonDegenerate());
    
    //////
    } 
    */      
  // have to slice if not zero
  //cerr << "convFunc nelements " << convFunc_p.nelements() << endl;
  if (convFunc_p.nelements() == 0) {

    Int npix = min(newAWConv.shape()[0], newAWConv.shape()[1]);
    //cerr << "####npix " << npix << " " << 2*max(awsupport)*oversamp_p << " oversamp " << oversamp_p << endl;
    if(npix < (2*max(awsupport+1)*oversamp_p)){
      npix=2*(max(awsupport)+1)*oversamp_p;
      //cerr << "aft npix " << npix << " shape " << newAWConv.shape() << endl;

      IPosition elshp=newAWConv.shape();
      elshp[0]=npix;
      elshp[1]=npix;
      convFunc_p=Array<Complex>(elshp,Complex(0.0));
      wgtConvFunc_p=Array<Complex>(elshp,Complex(0.0));
      MathUtils::putMiddle(convFunc_p, newAWConv);
      MathUtils::putMiddle(wgtConvFunc_p, newWtConv);
      
      
    }
    else{
      npix=2*(max(awsupport)+1)*oversamp_p;

      convFunc_p = MathUtils::getMiddle(newAWConv,  npix,  npix);
      wgtConvFunc_p = MathUtils::getMiddle(newWtConv,  npix,  npix);

    }
    convSizes_p.resize(trc[0] + 1);
    convSizes_p.set(npix);
    convSupport_p.resize(trc[0]+1);
    convSupport_p = awsupport.row(nfreqs-1);
  } else {
    // Appending PA changing only
    //in case support is bigger for this PA.

    IPosition newshp = convFunc_p.shape();

    if (newshp[0] < (2 * (max(awsupport) + 1) * oversamp_p)){
      //cerr << "@@@@RESHAPING " << endl;
      IPosition elshp = newshp;
      elshp[0] = (2 * (max(awsupport) + 1) * oversamp_p);
      elshp[1] = (2 * (max(awsupport) + 1) * oversamp_p);
      Array<Complex> tempC(elshp, Complex(0.0));
      Array<Complex> tempW(elshp,Complex(0.0));
      MathUtils::putMiddle(tempC, convFunc_p);
      MathUtils::putMiddle(tempW,wgtConvFunc_p);
      convFunc_p.reference(tempC);
      wgtConvFunc_p.reference(tempW);
    }
      // asumming same freqs for now
      newshp(4) = startrow + awConv.shape()(4);
    Int npix = min(newAWConv.shape()[0],  newAWConv.shape()[1]);
    if (npix <= newshp[0]) {
      IPosition blcadded(5, 0, 0, 0, 0, startrow);
      IPosition trcadded = newshp - 1;
      
      if(convFunc_p.shape()(4)< newshp[4]){
        convFunc_p.resize(newshp, True);  
        wgtConvFunc_p.resize(newshp, True);
      }
      convFunc_p(blcadded, trcadded).set(0.0);
      wgtConvFunc_p(blcadded, trcadded).set(0.0);
      Array<Complex> c = convFunc_p(blcadded, trcadded);
      MathUtils::putMiddle(c, newAWConv);
      Array<Complex> d=wgtConvFunc_p(blcadded, trcadded);
      MathUtils::putMiddle(d, newWtConv);
      
      convSizes_p.resize(trc[0]+1, true);
      convSizes_p(blc,  trc).set(newshp[0]);
      convSupport_p.resize(trc[0]+1, true);
      convSupport_p(blc, trc)= awsupport.row(nfreqs-1);
    } else {
      IPosition blcadded(5, 0, 0, 0, 0, startrow);
      if(newshp.product()>0 && newshp[0] < npix)
        npix=newshp[0];
      IPosition trcadded = newshp-1;
      if (convFunc_p.shape()(4) < newshp[4]) {
        convFunc_p.resize(newshp, True);
        wgtConvFunc_p.resize(newshp, True);
      }
      convFunc_p(blcadded, trcadded) =
          MathUtils::getMiddle(newAWConv, npix, npix);
      wgtConvFunc_p(blcadded, trcadded) =
            MathUtils::getMiddle(newWtConv, npix, npix);
      convSizes_p.resize(trc[0] + 1, true);
      convSizes_p(blc, trc).set(npix);
      convSupport_p.resize(trc[0] + 1, true);
      convSupport_p(blc, trc) = awsupport.row(nfreqs - 1);
    }
  }
  //cerr << "Shape at the end " << convFunc_p.shape() << endl;
}
Array<Complex>& AWConvFuncHolder::getConvFunc() {
  return convFunc_p;
  
}
Array<Complex>& AWConvFuncHolder::getWeightConvFunc() {
  return wgtConvFunc_p;
  
}
Vector<Int> AWConvFuncHolder::getConvSizes() {
  return convSizes_p;
}
Vector<Int> AWConvFuncHolder::getConvSupports() {
  
 return convSupport_p; 
}

Array<Complex> &AWConvFuncHolder::getConvFuncHPG() { return convFuncHPG_p; }
Array<Complex> &AWConvFuncHolder::getWeightConvFuncHPG() { return wgtConvFuncHPG_p; }
void AWConvFuncHolder::resetHPGConvFuncs(const vi::VisBuffer2 &vb){
  // A given set for hph will hold all w's
  wValsHPG_p.resize();
  wValsHPG_p = wVals_p;
  freqValsHPG_p.resize(freqVals_p.nelements(), False);
  Double fmin = min(vb.getFrequencies(0));
  Double fmax = max(vb.getFrequencies(0));
  uint indx = 0;
  vector<bool> usedfreq(freqVals_p.nelements());
  std::fill(usedfreq.begin(), usedfreq.end(), false);
  //cerr << "fmin " << fmin << " fmax " << fmax << "  freqs " << freqVals_p << "   " << (freqVals_p[0] >= fmin) << "   " <<(freqVals_p[0] <= fmax) << endl;
  if(freqVals_p.nelements() >1){
  for (uint k = 0; k < freqVals_p.nelements(); ++k) {
    if(freqVals_p[k] >= fmin && freqVals_p[k] <= fmax){
      freqValsHPG_p[indx] = freqVals_p[k];
      usedfreq[k] = true;
      ++indx;
    }
    
  }
  if(indx==0){ //some single channel spw will do this
    Double diffFreq=1e40;
    for (uint k = 0; k < freqVals_p.nelements(); ++k) {
      if(abs(fmax-freqVals_p[k]) < diffFreq){
        diffFreq=abs(fmax-freqVals_p[k]);
        freqValsHPG_p[0]=freqVals_p[k];
        std::fill(usedfreq.begin(),usedfreq.end(), false);
        usedfreq[k]=true;
        indx=1;
      }
    }

  }
  }
  else{
    //only one freq so it has to match
    usedfreq[0]=true;
    indx=1;

  }
  freqValsHPG_p.resize(indx, True);
  IPosition cshap = convFunc_p.shape();
  cshap[3] = indx;
  wgtConvFuncHPG_p.resize(cshap);
  convFuncHPG_p.resize(cshap);
  //cerr << "spw " << vb.spectralWindows()(0) << " CSHAP " << cshap << " orig " << convFunc_p.shape() << endl;
  IPosition blcin(5, 0);
  IPosition trcin = convFunc_p.shape() - 1;
  IPosition blcout(5,  0);
  IPosition trcout = cshap - 1;
  indx = 0;
  for (uint k = 0; k < freqVals_p.nelements(); ++k) {
    if(usedfreq[k]){
      blcin[3] = k;
      trcin[3] = k;
      blcout[3] = indx;
      trcout[3] = indx;
      for (uint j = 0; j < wVals_p.nelements(); ++j) {
        blcin[4] = j;
        trcin[4] = j;
        blcout[4] = j;
        trcout[4] = j;
        wgtConvFuncHPG_p(blcout, trcout) = wgtConvFunc_p(blcin, trcin);
        convFuncHPG_p(blcout, trcout) = convFunc_p(blcin, trcin);
      }
      ++indx;
    }
  }
}

/////////////////////
void AWConvFuncHolder::getConvFuncs(Vector<Int> &polMap, Vector<Int> &chanMap,
                                    Vector<Int> &rowMap,
                                    Array<Complex> &convFunc,
                                    Array<Complex> &wgtConvFunc,
                                    const vi::VisBuffer2 &vb,
                                    const Matrix<Double> &rotuvw,
                                    const Vector<Double> & interpFreqs,
                                    const Bool predictMode, 
                                    const bool ispsf) {

  Vector<Int> cmap;
  Vector<Int> pmap;
  Vector<Int> rmap;

  getConvIndices(pmap, cmap, rmap, vb, rotuvw, interpFreqs, predictMode, ispsf);
  //cerr << "pmap "<< pmap << endl;
  //cerr << "MIN Max rmap" << min(rmap) << "  " << max(rmap) << endl;
  std::vector<Int> pmapused = pmap.tovector();
  {
    std::sort(pmapused.begin(), pmapused.end());
    auto last = std::unique(pmapused.begin(), pmapused.end());
    pmapused.erase(last, pmapused.end());
  }
  std::vector<Int> cmapused = cmap.tovector();
  {
    std::sort(cmapused.begin(), cmapused.end());
    auto last = std::unique(cmapused.begin(), cmapused.end());
    cmapused.erase(last, cmapused.end());
  }
  std::vector<Int> rmapused = rmap.tovector();
  {
    std::sort(rmapused.begin(), rmapused.end());
    auto last = std::unique(rmapused.begin(), rmapused.end());
    rmapused.erase(last, rmapused.end());
  }
  {
    vector<Int> cpRmapUsed = rmapused;

    // lets move the -ve values to the end  -ve means -w which means we have to
    // conjugate the plane
    vector<int>::iterator it = remove_if(rmapused.begin(), rmapused.end(),
                                         [](const int i) { return i < 0; });
    rmapused.erase(it, rmapused.end());
    for (auto cit = cpRmapUsed.rbegin(); cit != cpRmapUsed.rend(); ++cit){
      if(*cit <0)
        rmapused.push_back(*cit);  
    }
  }
  //cerr << "#####rmapused " << rmapused << endl;
  IPosition shp(5, convFunc_p.shape()[0], convFunc_p.shape()[1],
                pmapused.size(), cmapused.size(), rmapused.size());
  polMap.resize(pmap.shape());
  for (uint j = 0; j < polMap.nelements(); ++j) {
    for (uint k = 0; k < pmapused.size(); ++k) {
      if (pmap[j] == pmapused[k]) {
        polMap[j] = k;
      }
    }
  }
  chanMap.resize(cmap.shape());
  //std::vector<int>cindex(cmapused.size());
  for (uint j = 0; j < chanMap.nelements(); ++j) {
    for (uint k = 0; k < cmapused.size(); ++k) {
      if (cmap[j] == cmapused[k]){
        chanMap[j] = k;
      }
    }
  }
  rowMap.resize(rmap.shape());
  for (uint j = 0; j < rowMap.nelements(); ++j) {
    for (uint k = 0; k < rmapused.size(); ++k) {
      if (rmap[j] == rmapused[k]){
        //rowmap is -ve for -ve w
        rowMap[j] = k;
      }
    }
  }
  //cerr << "old rmapused" << Vector<int>(rmapused) << " cmap " << Vector<int>(cmapused) << " pmap " << Vector<int>(pmapused) << endl;
  convFunc.resize(shp);
  wgtConvFunc.resize(shp);
  IPosition inblc(5, 0, 0, 0, 0, 0);
  IPosition intrc(5, shp[0] - 1, shp[1] - 1, 0, 0, 0);
  IPosition outblc(5, 0, 0, 0, 0, 0);
  IPosition outtrc(5, shp[0] - 1, shp[1] - 1, 0, 0, 0);

  for (uint r = 0; r < rmapused.size(); ++r){
    inblc[4] = abs(rmapused[r]);
    intrc[4] = abs(rmapused[r]);
    outblc[4] = r;
    outtrc[4] = r;
    for (uint c = 0; c < cmapused.size(); ++c) {
      inblc[3] = cmapused[c];
      intrc[3] = cmapused[c];
      outblc[3] = c;
      outtrc[3] = c;
      for (uint p = 0; p < pmapused.size(); ++p) {
        inblc[2] = pmapused[p];
        intrc[2] = pmapused[p];
        outblc[2] = p;
        outtrc[2] = p;
        // rowmap is -ve for -ve w
        convFunc(outblc, outtrc) = rmapused[r] > 0
                                       ? convFunc_p(inblc, intrc)
                                       : conj(convFunc_p(inblc, intrc));
        wgtConvFunc(outblc, outtrc) = wgtConvFunc_p(inblc, intrc);
      }
    }
  }
}
//////////////////////  
  

//////////////////////  

void AWConvFuncHolder::getConvIndices(Vector<Int>& polMap, Vector<Int>& chanMap, Vector<Int>& rowMap,  const vi::VisBuffer2& vb, const Matrix<Double>& rotuvw, const Vector<Double>& interpFreqs, 
  const Bool predictMode, const bool ispsf) {
  // Lets do the polmap
  Vector<Stokes::StokesTypes> visPolMap(vb.getCorrelationTypesSelected());
  polMap.resize(visPolMap.nelements());
  polMap.set(-1);
  if (!dosquint_p)
    polMap.set(0);
  else{
    for (uint k = 0; k < polMap.nelements(); ++k) {
     for (uint j = 0; j < polVals_p.nelements(); ++j) {
      if (visPolMap[k] == polVals_p[j])
        polMap[k] = j;
     }
    }
  }
  // Lets do chanMap
  chanMap.resize(interpFreqs.nelements());
  chanMap.set(-1);
  Vector<Double>visFreq = interpFreqs;
  for (uint k = 0; k < chanMap.nelements(); ++k) {
    Double minDiff = 1e40;
    Int indexF = -1;
    for (uint j = 0; j < freqVals_p.nelements(); ++j) {
      if (fabs(freqVals_p[j] -visFreq[k]) < minDiff) {
        minDiff = fabs(freqVals_p[j] -visFreq[k]);
        indexF = j;
      }
    }
    chanMap[k] = indexF;
        
  }
  // Now to the complicated rowMap
  Vector<Int> antPairIndex(vb.nRows(), -1);
  Vector<Int> wIndex(vb.nRows(), -1);
  Vector<Int> paIndex(vb.nRows(), -1);
  //Assuming pa is the same for this vb which is usually associated with time.
  // using utils.cc global function
  Double paval = refim::getPA(vb);
  if(predictMode){
    paval -= C::pi;
    paval = atan2(sin(paval), cos(paval));
  }
  Int tmpPAInd = -1;
  Double minDiff = 1e40;
  for (uint k = 0; k <paVals_p.nelements(); ++k) {
    if (fabs(paval-paVals_p[k]) < minDiff) {
        tmpPAInd = k;
        minDiff = fabs(paval-paVals_p[k]);
    }
  }
  paIndex.set(tmpPAInd);
  //cerr << "paVal " << paval << " predictMode " << predictMode << endl;
  // For antenna pairs ..for homogenous arrays only one pair is necessary
  if ( (antpairVals_p.nelements() == 1) && antpairVals_p[0] == std::pair<int,  int>(-1, -1)) {
      antPairIndex.set(0);
  }
  else{
      for (uint k = 0; k < vb.nRows(); ++k) {
       std::pair<int, int> antpair = std::make_pair(vb.antenna1()[k],  vb.antenna2()[k]);
       for (uint j = 0; j < antpairVals_p.nelements();++j ) {
        if (antpairVals_p[j] == antpair)
          antPairIndex[k] = j;
       }
         
      }
  }
  Double invlamda = mean(vb.getFrequencies(0))/C::c;
  for (uint k = 0; k < vb.nRows();++k) {
    minDiff = 1e40;
    Int tmpWInd = -1;
    Double w = ispsf ? 0 : rotuvw.row(2)[k] *invlamda;
    for (uint j =0; j < wVals_p.nelements();++j ) {
      if (fabs(fabs(w)-wVals_p[j]) < minDiff) {
       minDiff = fabs(fabs(w)-wVals_p[j]);
       tmpWInd = j;
      }
    }
    wIndex[k] = (w > 0)? -tmpWInd : tmpWInd;
  }
  // Now lets search for combination of all 3
  rowMap.resize(vb.nRows());
  for (uint k = 0; k < vb.nRows();++k) {
    for (uint j = 0; j < rowAxisWVals_p.nelements(); ++j) {
     if ( (abs(wIndex[k]) == rowAxisWVals_p[j]) && (paIndex[k] == rowAxisPAVals_p[j]) && (antPairIndex[k] == rowAxisAntennaPair_p[j]) ) {
      rowMap[k] = wIndex[k] >= 0 ? j : -j;
     }
    }
    
  }
  // A little dab will d'ya
  
}


void AWConvFuncHolder::getConvIndicesHPG(Vector<Int> &polMap, Vector<Int> &chanMap,
                                     Vector<Int> &rowMap,
                                     const vi::VisBuffer2 &vb,
                                     const Matrix<Double> &rotuvw) {
  Vector<Stokes::StokesTypes> visPolMap(vb.getCorrelationTypesSelected());
  polMap.resize(visPolMap.nelements());
  //HPG is doing I single plane gridding
  polMap.set(0);
  // Lets do chanMap 
  chanMap.resize(vb.nChannels());
  chanMap.set(-1);
  Vector<Double> visFreq = vb.getFrequencies(0);
  for (uint k = 0; k < chanMap.nelements(); ++k) {
    Double minDiff = 1e40;
    Int indexF = -1;
    for (uint j = 0; j < freqValsHPG_p.nelements(); ++j) {
      if (fabs(freqValsHPG_p[j] - visFreq[k]) < minDiff) {
        minDiff = fabs(freqValsHPG_p[j] - visFreq[k]);
        indexF = j;
      }
    }
    chanMap[k] = indexF;
  }
  //As there is no PA or antenna pair to deal with HPG...windex should be rowMap
  Vector<Int> wIndex(vb.nRows(), 0);
  Double invlamda = mean(vb.getFrequencies(0)) / C::c;
  for (uint k = 0; k < vb.nRows(); ++k) {
    Double minDiff = 1e40;
    Int tmpWInd = -1;
    Double w = rotuvw.row(2)[k] * invlamda;
    for (uint j = 0; j < wValsHPG_p.nelements(); ++j) {
      if (fabs(fabs(w) - wValsHPG_p[j]) < minDiff) {
        minDiff = fabs(fabs(w) - wValsHPG_p[j]);
        tmpWInd = j;
      }
    }
    wIndex[k] = tmpWInd;
  }
  rowMap.resize();
  rowMap = wIndex;
  //cerr << "FID " << vb.fieldId()(0) << " SPID " << vb.spectralWindows()(0) << " winDex " << wIndex << endl;
}

  Vector<Double> AWConvFuncHolder::getPointingPhaseShift(
      const vi::VisBuffer2 &vb, bool usePointingTable) {
    Bool hasValidPointing = False;
    if (vbutil_p.use_count() == 0)
      vbutil_p = std::make_shared<VisBufferUtil>(vb);
    MDirection ant1PointVal;
    if(Table::isReadable(vb.ms().pointingTableName())){
      hasValidPointing=usePointingTable &&  (vb.ms().pointing().nrow() >0);
    }
    DirectionCoordinate dc=outcsys_p.directionCoordinate(0);
   
    if(hasValidPointing){
      //ant1PointingCache_p[val]=vb.direction1()[0];
      ant1PointVal=vbutil_p->getPointingDir(vb, vb.antenna1()(0), 0, dc.directionType());
    }
    else
      ant1PointVal=vbutil_p->getPhaseCenter(vb);
    MSColumns mscol(vb.ms());
    String tel;
    if (vb.subtableColumns().observation().nrow() > 0) {
      tel =vb.subtableColumns().observation().telescopeName()(mscol.observationId()(0));
      }
    MEpoch::Types timeMType;
    casacore::Unit timeUnit;
    timeMType=MEpoch::castType(mscol.timeMeas()(0).getRef().getType());
    timeUnit=Unit(mscol.timeMeas().measDesc().getUnits()(0).getName());
    MPosition pos;
    MDirection dirOnImage;
    MeasTable::Observatory(pos,tel);
    //need to conver antpoint frame to image frame
    if(dc.directionType() !=  MDirection::castType(ant1PointVal.getRef().getType())){
    	
      MEpoch timenow(Quantity(vb.time()(0), timeUnit), timeMType);
      MeasFrame pointFrame(timenow, pos);
      MDirection::Ref elRef(dc.directionType(), pointFrame);
      dirOnImage=MDirection::Convert(ant1PointVal, elRef)();
      
    }
    else{
      dirOnImage=ant1PointVal;
    
    }
    
    Vector<Double> thePix(2);
    dc.toPixel(thePix, dirOnImage);
    //shift from center
    thePix(0) = thePix(0) - Double(nx_p / 2);
    thePix(1) = thePix(1) - Double(ny_p / 2);

    //phase gradient per pixel to apply
    thePix(0) = -thePix(0)*2.0*C::pi/Double(nx_p)/Double(oversamp_p);
    thePix(1) = -thePix(1) * 2.0 * C::pi / Double(ny_p) / Double(oversamp_p);
    //cerr << std::setprecision(12) << "fid " << vb.fieldId()(0) << " POINT shift " << thePix << endl;

    return thePix;

    
    
    
    
}

} // # namespace refim ends
}//namespace CASA ends

  
