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
    log_l << "Sum of weights: " << weights << " ";
    if (griddedData2.nelements() > 0) {
      log_l << max(griddedData2) << " " << min(griddedData2);
    };
    log_l << LogIO::POST;
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

    // cerr << "min correction " << min(sincConvX) << "    " << min(sincConvY)
    // << endl;
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
    Array<Complex> gwts;
    Bool removeDegenerateAxis = false;
    griddedWeights.get(gwts, removeDegenerateAxis);
    resampleCFToGrid(gwts, vbs, vb);
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
    Array<DComplex> gwts;
    Bool removeDegenerateAxis = false;
    griddedWeights_D.get(gwts, removeDegenerateAxis);
    resampleCFToGrid(gwts, vbs, vb);
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
  cerr << "###Max of imstore-> model " << max((imstore->model())->get())
       << endl;
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

}; // namespace refim
}; // namespace casa
