// -*- C++ -*-
//# AWProjectWBFTHPG.cc: Implementation of AWProjectWBFTHPG class
//# Copyright (C) 2021
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

#include <synthesis/TransformMachines2/AWProjectWBFTHPG.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <synthesis/ImagerObjects/SIImageStore.h>
#include <synthesis/ImagerObjects/SimpleSIImageStore.h>
#include <synthesis/TransformMachines/StokesImageUtil.h>
#include <synthesis/TransformMachines2/AWVisResampler.h>
#include <synthesis/TransformMachines2/SimplePBConvFunc.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/Slice.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/OS/HostInfo.h>
#include <casacore/casa/Utilities/CompositeNumber.h>
#include <casacore/images/Images/ImageInterface.h>
#include <casacore/images/Images/PagedImage.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <sstream>

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN
namespace refim {
//---------------------------------------------------------------
//
ImageInterface<Complex> &AWProjectWBFTHPG::getImage(Matrix<Float> &weights,
                                                    Bool normalize) {
  LogIO log_l(LogOrigin("AWProjectWBFTHPG", "getImage[R&D]"));
  //
  AlwaysAssert(image, AipsError);

  weights.resize(sumWeight.shape());
  convertArray(weights, sumWeight);
  //
  // If the weights are all zero then we cannot normalize otherwise
  // we don't care.
  //
  if (max(weights) == 0.0)
    log_l << "No useful data in " << name() << ".  Weights all zero"
          << LogIO::POST;
  else {
    /*log_l << "Sum of weights: " << weights << " ";
    if (griddedData2.nelements() > 0) {
      log_l << max(griddedData2) << " " << min(griddedData2);
    };
    log_l << LogIO::POST;*/
    cerr << "Sum of weights: " << setprecision(20) << weights << endl;
  }

  //
  // x and y transforms (lattice has the gridded vis.  Make the
  // dirty images)
  //
  if (useDoubleGrid_p) {
    // ArrayLattice<DComplex> darrayLattice(griddedData2);
    //  {
    //    griddedData.resize(griddedData2.shape());
    //    convertArray(griddedData, griddedData2);
    //    storeArrayAsImage(String("cgrid_"+visResampler_p->name()+".im"),
    //    image->coordinates(), griddedData);
    //  }
    // LatticeFFT::cfft2d(darrayLattice,false);

    griddedData.resize(griddedData2.shape());
    convertArray(griddedData, griddedData2);
    SynthesisUtilMethods::getResource("mem peak in getImage");

    // Don't need the double-prec grid anymore...
    griddedData2.resize();
    lattice = new ArrayLattice<Complex>(griddedData);
  } else {
    lattice = new ArrayLattice<Complex>(griddedData);
    // // cerr << "##### " << griddedData2.shape() << endl;
    // lattice=arrayLattice;
    // LatticeFFT::cfft2d(*lattice,false);
  }

  //
  // Now normalize the dirty image.
  //
  // Since *lattice is not copied to *image till the end of this
  // method, normalizeImage also needs to work with Lattices
  // (rather than ImageInterface).
  //
  // //normalizeImage(*lattice,sumWeight,*avgPB_p,fftNormalization);
  //	normalizeImage(*lattice,sumWeight,*avgPB_p, *avgPBSq_p,
  //fftNormalization);

  // nx ny normalization from GridFT...
  {
    Int inx = lattice->shape()(0);
    Int iny = lattice->shape()(1);
    Vector<Complex> correction(inx);
    correction = Complex(1.0, 0.0);
    Vector<Float> sincConvX(inx);
    for (Int ix = 0; ix < inx; ix++) {
      Float x = C::pi * Float(ix - inx / 2) / (Float(nx) * Float(convSampling));
      if (ix == inx / 2) {
        sincConvX(ix) = 1.0;
      } else {
        sincConvX(ix) = sin(x) / x;
      }
    }
    Vector<Float> sincConvY(iny);
    for (Int ix = 0; ix < iny; ix++) {
      Float x = C::pi * Float(ix - iny / 2) / (Float(ny) * Float(convSampling));
      if (ix == iny / 2) {
        sincConvY(ix) = 1.0;
      } else {
        sincConvY(ix) = sin(x) / x;
      }
    }

     //cerr << "NORM " << normalize << " min correction " << min(sincConvX) << "    " << min(sincConvY) << endl;

    //  Do the Grid-correction
    IPosition cursorShape(4, inx, 1, 1, 1);
    IPosition axisPath(4, 0, 1, 2, 3);
    LatticeStepper lsx(lattice->shape(), cursorShape, axisPath);
    LatticeIterator<Complex> lix(*lattice, lsx);

    for (lix.reset(); !lix.atEnd(); lix++) {
      Int pol = lix.position()(2);
      Int chan = lix.position()(3);
      if (weights(pol, chan) != 0.0) {
        Int iy = lix.position()(1);
        for (Int ix = 0; ix < nx; ix++)
          correction(ix) = 1 / (sincConvX(ix) * sincConvY(iy));
        // cerr <<"Min max correction " << min(correction) << "     " <<
        // max(correction) << endl;
        lix.rwVectorCursor() *= correction;
        if (normalize) {
          Complex rnorm(Float(inx) * Float(iny) / weights(pol, chan));
          lix.rwCursor() *= rnorm;
        }
      } else {
        lix.woCursor() = 0.0;
      }
    }

    // for(lix.reset();!lix.atEnd();lix++)
    //   {
    //     Int pol=lix.position()(2);
    //     Int chan=lix.position()(3);
    //     if (normalize)
    // 	{
    // 	  if(weights(pol,chan)!=0.0)
    // 	    {
    // 	      Complex rnorm(Float(inx)*Float(iny)/(sincConv(inx)*sincConv(iny)*
    // weights(pol,chan) )); 	      lix.rwCursor()*=rnorm;
    // 	    }
    // 	  else
    // 	    lix.woCursor()=0.0;
    // 	}
    //     else
    // 	lix.rwCursor() /= sincConv(inx)*sincConv(iny);
    //   }
  }
  if (!isTiled) {
    //
    // Check the section from the image BEFORE converting to a lattice
    //
    IPosition blc(4, (nx - image->shape()(0) + (nx % 2 == 0)) / 2,
                  (ny - image->shape()(1) + (ny % 2 == 0)) / 2, 0, 0);
    IPosition stride(4, 1);
    IPosition trc(blc + image->shape() - stride);
    //
    // Do the copy
    //
    //cerr << "blc" << blc << " trc " << trc << " min max " << min(griddedData) << "  max " << max(griddedData) << endl;
    image->put(griddedData(blc, trc));

    if (!lattice.null())
      lattice = 0;
    griddedData.resize(IPosition(1, 0));
  }
  {
    // TempImage<Complex> tt(lattice->shape(), image->coordinates());
    // tt.put(lattice->get());
    // storeImg(String("uvgrid"+visResampler_p->name()+".im"), *image,true);
  }

  return *image;
}
//////////////////////////////////////////

void AWProjectWBFTHPG::getWeightImage(ImageInterface<Float> &weightImage,
                                      Matrix<Float> &weights) {
  /// This is a mess nothing is initialized properly..
  // lets make a guess :)
  // SB: Review the guess below (not by SB) to see if it is still
  // required after the code cleanup.
  if (avgPB_p && (avgPB_p->shape()).size() == 0)
    avgPB_p = nullptr;
  if (!avgPB_p) {
    (getImage(weights, false));
    IPosition cursorShape(4, image->shape()[0], image->shape()[1], 1, 1);
    IPosition axisPath(4, 0, 1, 2, 3);
    LatticeStepper lsx(image->shape(), cursorShape, axisPath);
    LatticeIterator<Complex> lix(*image, lsx);
    for (lix.reset(); !lix.atEnd(); lix++) {
      Int pol = lix.position()(2);
      Int chan = lix.position()(3);
      if (weights(pol, chan) != 0.0) {
        Complex rnorm(1 / weights(pol, chan));
        lix.rwCursor() *= rnorm;
      } else {
        lix.rwCursor() = 0;
      }
    }
    StokesImageUtil::ToStokesPSF(weightImage, *image);
    setWeightImage(weightImage);
  } else {
    weightImage.resize(avgPB_p->shape());
    weightImage.copyData(*avgPB_p);
    weights.resize(sumWeight.shape());
    convertArray(weights, sumWeight);
  }
  avgPBReady_p = True;
}
//
//---------------------------------------------------------------
// Methods to accumulate data on the grid.  These trigger gridding
// to make *one* type of image: weights, PSF or residual depending
// on the vbs.ftmType_p. This is required with HPG since it grids
// for only one type of image at a time for efficiency.  This is
// different from the pattern on the CPU where weights are gridded
// along with the gridding for the first image (PSF or residual).
// Hence the specialization here.
//
void AWProjectWBFTHPG::resampleDataToGrid(Array<Complex> &griddedData_l,
                                          VBStore &vbs, const VisBuffer2 &vb,
                                          Bool &dopsf) {
  if (vbs.ftmType_p != casa::refim::FTMachine::WEIGHT)
    AWProjectFT::resampleDataToGrid(griddedData_l, vbs, vb, dopsf);
  // if (!avgPBReady_p)
  else {
    //
    // Get a reference to the pixels of griddedWeights (a
    // TempImage!)
    //
    //Array<Complex> gwts;
    //Bool removeDegenerateAxis = false;
    //griddedWeights.get(gwts, removeDegenerateAxis);
    //resampleCFToGrid(gwts, vbs, vb);
    vbs.ftmType_p=casa::refim::FTMachine::WEIGHT;  
    Int nDataChan = vbs.flagCube_p.shape()[1];
    vbs.startChan_p = 0; vbs.endChan_p = nDataChan;
    Bool locdopsf=true;
     AWProjectFT::resampleDataToGrid(griddedData_l, vbs, vb, locdopsf);
    
  }
};
//
//---------------------------------------------------------------
//
void AWProjectWBFTHPG::resampleDataToGrid(Array<DComplex> &griddedData_l,
                                          VBStore &vbs, const VisBuffer2 &vb,
                                          Bool &dopsf) {
  if (vbs.ftmType_p != casa::refim::FTMachine::WEIGHT)
    AWProjectFT::resampleDataToGrid(griddedData_l, vbs, vb, dopsf);
  // if (!avgPBReady_p)
  else {
    //
    // Get a reference to the pixels of griddedWeights (a
    // TempImage!)
    //
    
    //Array<DComplex> gwts;
    //Bool removeDegenerateAxis = false;
    //griddedWeights_D.get(gwts, removeDegenerateAxis);
    //resampleCFToGrid(griddedData_l, vbs, vb);
    vbs.ftmType_p=casa::refim::FTMachine::WEIGHT;  
    Int nDataChan = vbs.flagCube_p.shape()[1];
    vbs.startChan_p = 0; vbs.endChan_p = nDataChan;
    Bool locdopsf=true;
    AWProjectFT::resampleDataToGrid(griddedData_l, vbs, vb, locdopsf);
  }
};
//
//---------------------------------------------------------------
//
void AWProjectWBFTHPG::initializeToVisNew(const VisBuffer2 &vb,
                                          CountedPtr<SIImageStore> imstore) {

  Matrix<Float> tempWts;

  if (!(imstore->forwardGrid()).get())
    throw(AipsError("FTMAchine::InitializeToVisNew error imagestore has no "
                    "valid grid initialized"));
  // Convert from Stokes planes to Correlation planes
  LatticeLocker lock1(*(imstore->model()), FileLocker::Read);
  // cerr << "###Max of imstore-> model " << max((imstore->model())->get())
  //      << endl;
  stokesToCorrelation(*(imstore->model()), *(imstore->forwardGrid()));

  if (vb.polarizationFrame() == MSIter::Linear) {
    StokesImageUtil::changeCStokesRep(*(imstore->forwardGrid()),
                                      StokesImageUtil::LINEAR);
  } else {
    StokesImageUtil::changeCStokesRep(*(imstore->forwardGrid()),
                                      StokesImageUtil::CIRCULAR);
  }
  setFTMType(refim::FTMachine::RESIDUAL);
  visResampler_p->setModelImage((imstore->forwardGrid()));
}
//-------------------------------------------------------------------------
  //  
  void AWProjectWBFTHPG::setupVBStore(VBStore& vbs,
				 const VisBuffer2& vb, 
				 const Matrix<Float>& imagingweight,
				 const Cube<Complex>& visData,
				 const Matrix<Double>& uvw,
				 const Cube<Int>& flagCube,
				 const Vector<Double>& dphase,
				 const Bool& dopsf,
				 const Vector<Int>& /*gridShape*/)
  {
    vbs.vb_p = &vb;
    vbs.wbAWP_p=wbAWP_p;
    vbs.ftmType_p=ftmType_p;
    vbs.nWPlanes_p = nWPlanes_p;
    //cerr << "HPG setupvbstore " << endl;

    visResampler_p->setParams(uvScale,uvOffset,dphase);
    visResampler_p->setMaps(chanMap, polMap);
    
    //
    // Set up VBStore object to point to the relavent info. of the VB.
    //
    vbs.imRefFreq_p = imRefFreq_p;
    vbs.nRow_p = vb.nRows();
    vbs.beginRow_p = 0;
    vbs.endRow_p = vbs.nRow_p;
    vbs.spwID_p = vb.spectralWindows()(0);
    vbs.nDataPol_p  = flagCube.shape()[0];
    vbs.nDataChan_p = flagCube.shape()[1];

    vbs.antenna1_p.reference(vb.antenna1());
    vbs.antenna2_p.reference(vb.antenna2());
    //vbs.paQuant_p = Quantity(getPA(vb),"rad");

    vbs.corrType_p.reference(vb.correlationTypes());

    vbs.uvw_p=uvw;
    vbs.imagingWeight_p.reference(imagingweight);
    vbs.visCube_p.reference(visData);

    vbs.freq_p.reference(vb.getFrequencies(0));

    vbs.rowFlag_p.reference(vb.flagRow());
    if(!usezero_p) 
      for (Int rownr=0; rownr<vbs.nRow_p; rownr++) 
	if(vb.antenna1()(rownr)==vb.antenna2()(rownr)) vbs.rowFlag_p(rownr)=true;

    vbs.flagCube_p.resize(flagCube.shape());  vbs.flagCube_p = false; vbs.flagCube_p(flagCube!=0) = true;
      
    vbs.conjBeams_p=conjBeams_p;

    
    // The following code is required only for GPU or multi-threaded
    //gridder.  Currently does not work without the rest of the
    //GPU/multi-threaded infrastructure (though, I (SB) thought this
    //was designed to be benign for normal gridding).
    //
    //visResampler_p->initializeDataBuffers(vbs);
  }
   void AWProjectWBFTHPG::init(const vi::VisBuffer2& vb) 
  {
    LogIO log_l(LogOrigin("AWProjectFT2", "init[R&D]"));

    nx    = image->shape()(0);
    ny    = image->shape()(1);
    npol  = image->shape()(2);
    nchan = image->shape()(3);
    
    
    sumWeight.resize(npol, nchan);
    sumCFWeight.resize(npol, nchan);
    
    wConvSize=max(1, nWPlanes_p);
    
    CoordinateSystem cs=image->coordinates();
    uvScale.resize(3);
    uvScale=0.0;
    uvScale(0)=Float(nx)*cs.increment()(0); 
    uvScale(1)=Float(ny)*cs.increment()(1); 
    uvScale(2)=Float(wConvSize)*abs(cs.increment()(0));
    
    Int index= cs.findCoordinate(Coordinate::SPECTRAL);
    SpectralCoordinate spCS = cs.spectralCoordinate(index);
    imRefFreq_p = spCS.referenceValue()(0);
    double f1, f2;
    spCS.toWorld(f1, double(-0.5));
    spCS.toWorld(f2, double(nchan)-0.5);
    auto frange=std::make_pair(f1, f2);
    uvOffset.resize(3);
    uvOffset(0)=nx/2;
    uvOffset(1)=ny/2;
    uvOffset(2)=0;
    
    if(gridder) delete gridder;
    gridder=0;
    gridder = new ConvolveGridder<Double, Complex>(IPosition(2, nx, ny),
						   uvScale, uvOffset,
						   "SF");
    makingPSF = false;
    


///We'll always use oversampling of 4 for GPU gridder
    //convSampling=4;
  //TESTOO
   
    convSampling=10;
    if (min(nx, ny) > 200)
      convSampling = 4;
    ///////

    // we are not doing parallactiv angle here
    Double painc=2*C::pi;

  if(awConvs_p.use_count()==0){
     String observatory=(vb.subtableColumns().observation()).telescopeName()(0);
    awConvs_p=std::make_shared<AWConvFuncHolder>((*image).coordinates(), nx, ny, 
                   False, painc, observatory, convSampling);
    vi::VisibilityIterator2 *vi= const_cast<VisibilityIterator2 *>(vb.getVi());
    
    
    std::vector<Double> freqs;
    std::set<Int> fields;
    std::vector<Double> pAs={0.0};
    //int validspw=-1;
    Double maxW=0.0;
    for (vi->originChunks(); vi->moreChunks(); vi->nextChunk()) {
          for (vi->origin(); vi->more(); vi->next()) {
              std::vector<Double> chunkfreq;
              fields.insert(vb.fieldId()(0));
              //matchChannel(vb);
              //cerr << "MAX chanMap" << chanMap << endl;
              //if (max(chanMap) > -1) 
              {

                SimplePBConvFunc::findUsefulChannels(chunkfreq, vb, frange);
                //cerr << "vbnchan " << vb.nChannels() << "chunkfreq "
                //     << chunkfreq << endl;
                if (chunkfreq.size() > 0) {
                  // validspw=vb.spectralWindows()(0);
                  // cerr << "SPW " << vb.spectralWindows()(0) << " freqs " <<
                  // Vector<Double>(chunkfreq) << endl;
                  std::move(chunkfreq.begin(), chunkfreq.end(),
                            std::back_inserter(freqs));
                  double maxfreqused =
                      *(std::max_element(chunkfreq.begin(), chunkfreq.end()));
                  if (nWPlanes_p > 1) {
                    // 	maxW=max(maxW,
                    // max(abs(vb.uvw().row(2)*max(vb.getFrequencies(0))))/C::c);
                    maxW = max(maxW,
                               max(abs(vb.uvw().row(2) * maxfreqused)) / C::c);
                  }
                }
              }
          }
    }
  
    
    //return vi to origin
    vi->originChunks(); vi->origin();
    
    std::sort(freqs.begin(),  freqs.end());
    auto last = std::unique(freqs.begin(),  freqs.end());
    freqs.erase(last,  freqs.end());
    // tell holder it is a single field or not
    (*awConvs_p).setSingleField((fields.size() == 1));

    if (nWPlanes_p == 0)
      nWPlanes_p = 1;
    Vector<Double> wVals(nWPlanes_p,0);
    if(nWPlanes_p >1){
      Double st=maxW/(Double(nWPlanes_p-1)*Double(nWPlanes_p-1));
      for (int k=0; k <nWPlanes_p; ++k)
        wVals[k]=Double(k*k)*st;
    }
    //cerr << "XXXXXINIT wVals " << wVals << endl;
    //cerr << "XXXXXfreqs " << Vector<Double>(freqs) << endl; 
    (*awConvs_p).addConvFunc(Vector<Double>(freqs), wVals, 0.0);
    /////TESTOO
   /*{
		Vector<Double>pixW(wVals.nelements());
		indgen(pixW);
		CoordinateSystem fiveAxis=image->coordinates();
    Vector<Int> stoks(4);
    stoks(0) = Stokes::RR;
    stoks(1) = Stokes::RL;
    stoks(2) = Stokes::LR;
    stoks(3) = Stokes::LL;
    StokesCoordinate stokesCoords(stoks);
    fiveAxis.replaceCoordinate(stokesCoords, 1);

		TabularCoordinate tab(pixW, wVals, "m", "W");
		fiveAxis.addCoordinate(tab);
		PagedImage<Complex> noo((awConvs_p->getConvFunc()).shape(), fiveAxis, "AWConvVals_"+String::toString(validspw));
		noo.put((awConvs_p->getConvFunc()));
	
	}*/


    ///////////////////////////
    visResampler_p->setConvFunc(awConvs_p);
  }
}
 void AWProjectWBFTHPG::initializeToSky(ImageInterface<Complex>& iimage,
				   Matrix<Float>& weight,
				   const VisBuffer2& vb)
  {
    LogIO log_l(LogOrigin("AWProjectWBFT2","initializeToSky[R&D]"));
  //  Timer tim;
  //  tim.mark();
    
    AWProjectFT::initializeToSky(iimage,weight,vb);
  //  cerr << "$$$$$ AWP initializesky " << tim.real()<< endl;

  //  tim.mark();
  init(vb);
  //cerr  << "$$$$$$$$$$$ init(vb) " << tim.real() << endl;

  }
  void AWProjectWBFTHPG::findConvFunction(const ImageInterface<Complex>& image,
				     const VisBuffer2& vb)
  {
 
 
  //  CoordinateSystem ftcoords;
  //We have to make sure awh is loaded 
    if(awConvs_p.use_count()==0)
      throw(AipsError("Programmer's error:Convolution function has not been set"));
/* 
    if(avgPBReady_p){
        LatticeExprNode le( max( *avgPB_p ) );
        Float avgPB_max=le.getFloat();
        
        if(avgPB_max <= 0.0) avgPBReady_p = false;
    }
    
    if(!avgPBReady_p) makeSensitivityImage(vb,image,*avgPB_p);

	
    
    verifyShapes(avgPB_p->shape(), image.shape());
*/
    
	
      
  }
}; // namespace refim
}; // namespace casa
