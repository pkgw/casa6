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
#include <iomanip>

namespace casa {  //# CASA namespace
namespace refim { //# namespace refactor imaging

using namespace casacore;
using namespace casa;
using namespace casa::refim;
using namespace std;

AWConvFuncHolder::AWConvFuncHolder(const CoordinateSystem& csys, const int nx, const int ny, const bool dosquint, const double paInc,  const String& obs,  const int& oversamp): painc_p(paInc),  dosquint_p(dosquint), outcsys_p(csys),  nx_p(nx),  ny_p(ny),  oversamp_p(oversamp) {
  
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
    convSizes_p = other.convSizes_p;
    convSupport_p.resize();
    convSupport_p = other.convSupport_p;
    aterm_p = other.aterm_p;
    
    
    
    
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
  //cerr << "@@@@@@@@@@@@@@@@@@@@@@wVals to be calc " << wVals_p << endl;
 } else if (wVals_p.nelements() != wVals.nelements())
   throw(AipsError("Cannot change W terms length right now"));
 //cerr << "########WVals_p " << wVals_p << endl;
 if (!dosquint_p) {
   paVals_p.resize(1);
   paVals_p[0] = 0.0;
 } else {
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
  Array<Complex> aWConv;
  Array<Complex> aWwtconv;
  Matrix<Int> awSupport;
  calcCsys_p = outcsys_p;
  calcNpix_p = min(nx_p,  ny_p);
  //cerr << "PAVALS " <<  paVals_p <<  " dosquint " << dosquint_p <<  endl;
  //cerr << "FREQS " << freqsToCalc << endl;
  for (uint k=0; k<paVals_p.nelements(); ++k){
    a.makeAWConvFunc(aWConv, aWwtconv,calcCsys_p,awSupport, calcNpix_p, freqsToCalc, wVals_p, dosquint_p, paVals_p[k]);
    cerr << "######MAX awsupp " << max(awSupport) << endl;
                                                   
    //append arrays and indices  
    appendConvFuncs(aWConv,  aWwtconv,  awSupport,  freqsToCalc,  paVals_p[k]);
  }
  
  
 
 
 return true;;
}
void AWConvFuncHolder::appendConvFuncs(const Array<Complex>& awConv,  const Array<Complex>& aWwtConv,  const Matrix<Int>& awsupport,  const Vector<Double>& newFreqs,  const Double paVal) {
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
  rowAxisAntennaPair_p.resize(trc[0]+1,  True);
  rowAxisAntennaPair_p(blc, trc).set(0);
  /// Let us rescale convolution function to match nx, ny and incr of image that makes 
  /// uvgrid
  Float factorX=fabs(calcCsys_p.increment()(0)/outcsys_p.increment()(0));
  Float factorY=fabs(calcCsys_p.increment()(1)/outcsys_p.increment()(1));
//  cerr <<  "####Factor " <<  factorX <<  "   " <<  factorY <<  endl;
  factorX = Float(nx_p) *Float(oversamp_p)/Float(calcNpix_p)/factorX;
  factorY = Float(ny_p) *Float(oversamp_p)/Float(calcNpix_p)/factorY;
//  cerr <<  "factors " <<  factorX <<  "   " <<  factorY <<  "nx,  ny" <<  nx_p << "   " << ny_p << " calcNpix " << calcNpix_p << " oversamp " << oversamp_p << endl;
  MathUtils m;
  Array<Complex>newAWConv = m.resampleViaFFT(awConv,  factorX,  factorY);
  Array<Complex> newWtConv = m.resampleViaFFT(aWwtConv,  factorX,  factorY);
  Float correcfac = float(awConv.shape()(0) *awConv.shape()(1) *oversamp_p*oversamp_p)/float(newAWConv.shape()(0) *newAWConv.shape()(1));
  //cerr <<  "correcfac " <<  correcfac  <<  "  "  <<  1.0/correcfac  <<  endl;
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
    }*/       
  // have to slice if not zero
  if (convFunc_p.nelements() == 0) {
    Int npix = min(newAWConv.shape()[0],  newAWConv.shape()[1]);
    //cerr << "npix " << npix << " " << 2*max(awsupport)*oversamp_p << endl;
    if(npix <= 2*max(awsupport)*oversamp_p){
      npix=2*(max(awsupport)+1)*oversamp_p;
      cerr << "aft npix " << npix << endl;
      IPosition elshp=newAWConv.shape();
      elshp[0]=npix;
      elshp[1]=npix;
      convFunc_p=Array<Complex>(elshp,Complex(0.0));
      wgtConvFunc_p=Array<Complex>(elshp,Complex(0.0));
      MathUtils::putMiddle(convFunc_p, newAWConv);
      MathUtils::putMiddle(wgtConvFunc_p, newWtConv);
      
      
    }
    else{
      convFunc_p = MathUtils::getMiddle(newAWConv,  npix,  npix);
      wgtConvFunc_p = MathUtils::getMiddle(newWtConv,  npix,  npix);
    }
    convSizes_p.resize(trc[0]+1);
    convSizes_p.set(npix);
    convSupport_p.resize(trc[0]+1);
    convSupport_p = awsupport.row(nfreqs-1);
  }
  else{
    // Appending PA changing only
    IPosition newshp = convFunc_p.shape();
    // asumming same freqs for now
    newshp(4) = newshp(4)+awConv.shape()(4);
    Int npix = min(newAWConv.shape()[0],  newAWConv.shape()[1]);
    if (npix !=  newshp[0])
    {
      
     cerr <<  "npix is not the same for a different PA" <<  endl; 
    }
    else{
      IPosition blcadded(5,  0,  0,  0,  0, convFunc_p.shape()[4]);
      IPosition trcadded = newshp-1;
      convFunc_p.resize(newshp,  True);
      wgtConvFunc_p.resize(newshp,  True);
      convFunc_p(blcadded,  trcadded) = MathUtils::getMiddle(newAWConv,  npix,  npix);
      wgtConvFunc_p(blcadded, trcadded) = MathUtils::getMiddle(newWtConv,  npix,  npix);
      convSizes_p.resize(trc[0]+1);
      convSizes_p(blc,  trc).set(npix);
      convSupport_p.resize(trc[0]+1);
      convSupport_p(blc, trc)= awsupport.row(nfreqs-1);
      
      
    }
  }
  
  
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

/////////////////////
void AWConvFuncHolder::getConvFuncs(Vector<Int> &polMap,Vector<Int> &chanMap,Vector<Int> &rowMap, Array<Complex> &convFunc, Array<Complex> &wgtConvFunc, 
                  const vi::VisBuffer2 &vb,
                  const Matrix<Double> &rotuvw){

  Vector<Int> cmap;
  Vector<Int> pmap;
  Vector<Int> rmap;
  getConvIndices(pmap, cmap, rmap, vb, rotuvw);
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
    vector<Int> cpRmapUsed=rmapused;
  //lets move the -ve values to the end  -ve means -w which means we have to conjugate the plane
   vector<int>::iterator it =
      remove_if(rmapused.begin(), rmapused.end(), [](const int i) { return i < 0; });
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
      if (pmap[j]==pmapused[k]){
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
        //rowmap is -ve for -ve w
        convFunc(outblc, outtrc) = rmapused[r] >0 ? convFunc_p(inblc, intrc) : conj(convFunc_p(inblc, intrc));
        wgtConvFunc(outblc, outtrc) = wgtConvFunc_p(inblc, intrc);
      }
    }
  }
}
//////////////////////  
  

void AWConvFuncHolder::getConvIndices(Vector<Int>& polMap, Vector<Int>& chanMap, Vector<Int>& rowMap,  const vi::VisBuffer2& vb, const Matrix<Double>& rotuvw) {
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
  chanMap.resize(vb.nChannels());
  chanMap.set(-1);
  Vector<Double>visFreq = vb.getFrequencies(0);
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
  Int tmpPAInd = -1;
  Double minDiff = 1e40;
  for (uint k = 0; k <paVals_p.nelements(); ++k) {
    if (fabs(paval-paVals_p[k]) < minDiff) {
        tmpPAInd = k;
        minDiff = fabs(paval-paVals_p[k]);
    }
  }
  paIndex.set(tmpPAInd);
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
    Double w = rotuvw.row(2)[k] *invlamda;
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
      rowMap[k] = wIndex[k] > 0 ? j : -j;
     }
    }
    
  }
  // A little dab will d'ya
  
}
} // # namespace refim ends
}//namespace CASA ends

  
