//# tALMAAperture.cc: This program tests the ALMAAperture class
//# Copyright (C) 1998,1999,2000
//# Associated Universities, Inc. Washington DC, USA.
//# Copyright (C) 2011 ESO (in the framework of the ALMA collaboration)
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$

//# Includes
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/aips.h>
#include <casacore/coordinates/Coordinates.h>
#include <casacore/measures/Measures.h>
#include <casacore/measures/Measures/MDirection.h>

#include <casacore/images/Images/PagedImage.h>
#include <casacore/images/Images/SubImage.h>
#include <casacore/images/Regions/ImageRegion.h>

#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Arrays/Slicer.h>
#include <casacore/lattices/LRegions/LCRegion.h>
#include <casacore/lattices/LRegions/LCSlicer.h>
#include <casacore/lattices/Lattices/LatticeIterator.h>
#include <casacore/lattices/Lattices/LatticeStepper.h>
#include <casacore/lattices/Lattices/TiledShape.h>

#include <components/ComponentModels/ComponentType.h>
#include <components/ComponentModels/Flux.h>
#include <components/ComponentModels/PointShape.h>
#include <components/ComponentModels/SkyComponent.h>
#include <components/ComponentModels/SpectralIndex.h>

#include <casacore/casa/BasicSL/String.h>
#include <synthesis/TransformMachines/PBMath.h>
#include <synthesis/TransformMachines/SkyJones.h>

#include <msvis/MSVis/VisBuffer.h>
#include <msvis/MSVis/VisSet.h>
#include <msvis/MSVis/VisibilityIterator.h>

#include <synthesis/Utilities/FFT2D.h>
#include <synthesis/TransformMachines2/EVLAAperture.h>
#include <synthesis/TransformMachines2/WPConvFunc.h>
#include <synthesis/TransformMachines2/AWConvFunc.h>
#include <synthesis/TransformMachines2/Utils.h>
#include <synthesis/TransformMachines/PBMath1DAiry.h>
#include <casacore/casa/OS/Timer.h>

using namespace std;
using namespace casa;
using namespace casacore;

casa::refim::MathUtils mUtils;

int findConvSize(CoordinateSystem& cs){
  Int index = cs.findCoordinate(Coordinate::SPECTRAL);
  Double freqVal = 0.0;
  cs.spectralCoordinate(index).toWorld(freqVal, 0.0);
  Double beamSize=C::c/freqVal/25.0;
  index = cs.findCoordinate(Coordinate::DIRECTION);
  Quantity incr(fabs(cs.directionCoordinate(index).increment()[0]), cs.directionCoordinate(index).worldAxisUnits()[0]);
  cerr << beamSize << " / " << incr << " = " <<   (beamSize/incr.getValue("rad")) << endl;
  int size=16*std::round((beamSize/(incr.getValue("rad")))/8);
  return size;
}
bool fakeConv(ImageInterface<Complex>& conv, ImageInterface<Complex>& out, Float factor, Int origSupp){
  
     IPosition convshp=conv.shape();
     IPosition outshp=out.shape();
     IPosition blc(4, convshp[0]/2-origSupp+1,  convshp[1]/2-origSupp+1, 0, 0);
     IPosition shp(4, 2*origSupp, 2*origSupp, convshp[2], convshp[3]);
     Slicer sl(blc, shp);
     Array<Complex> origConv=conv.getSlice(sl);
     Array<Complex> newConv;
     newConv = mUtils.resample(origConv, factor, factor);
     shp=newConv.shape();
     cerr << "newconv shape " << newConv.shape() << endl;
     blc=IPosition(4, outshp[0]/2-shp[0]/2+1,  outshp[1]/2-shp[1]/2+1, 0, 0);
     
     sl=Slicer(blc, shp);
     out.set(0.0);
     out.putSlice(newConv, blc);
     FFT2D ftim;
      ftim.c2cFFT(out, False);
     
      out.copyData(LatticeExpr<Complex> (out/max(abs(out))));
     return true;
  
}

bool resampleCopy(ImageInterface<Complex>& outImage, ImageInterface<Complex>& inImage, const Float factor){
  Array<Complex> inArr=inImage.get();
  Array<Complex> outArr;
  outArr = mUtils.resample(inArr, factor, factor);
  IPosition blcin(4, 0, 0, 0, 0);
  IPosition shpin=outArr.shape();
  IPosition blcout(4, 0, 0, 0, 0);
  for(int j=0; j <2; ++j){
    if(outArr.shape()[j] > outImage.shape()[j]){
      blcin[j]= (outArr.shape()[j]-outImage.shape()[j])/2;
      shpin[j]=outImage.shape()[j];
    }
    else if (outImage.shape()[j] > outArr.shape()[j]){
      blcout[j]= (outImage.shape()[j]-outArr.shape()[j])/2;
    }
  }
  Slicer insl(blcin, shpin);
  cerr << "outArr.shape " << outArr.shape() << " blcin " <<  blcin << " shpin " << shpin << " blcout " << blcout << endl; 
  outImage.putSlice(outArr(insl), blcout);
  
 return true; 
}

int main() {
  try {
    /// BYPASSING test for now
    return 0;
    // Form coordinate systems and test images
    CoordinateSystem coordsys;
    CoordinateSystem coordsys3;
    CoordinateSystem coordsys3Big;
    CoordinateSystem coordsys3Small;
    CoordinateSystem coordsys4;
    Double leFreq=2.029e9;
    cout << "new EVLAAperture" << endl;
    refim::EVLAAperture aa;
    aa.cacheVBInfo("EVLA", 25.0);
    cout << aa.name() << endl;
    Int bandid=aa.getBandID(leFreq);
    String bandname=aa.getBandName();
    Quantity fov(0.024,"rad");
    if(bandname=="EVLA_S")
      fov /=2;
    Quantity cell=fov/256;
    Quantity imcell(0.6, "arcsec");
    {
      Matrix<Double> xform(2, 2);
      xform = 0.0;
      xform.diagonal() = 1.0;
      Quantity RA; Quantity::read(RA, "13:03:44.242350");
      Quantity DEC; Quantity::read(DEC,"-10.55.12.00000");
      DirectionCoordinate dirCoords(
          MDirection::J2000, Projection(Projection::SIN), RA.get("rad").getValue(),
          DEC.get("rad").getValue(), imcell.get("rad").getValue(),
          imcell.get("rad").getValue(), xform, 2600, 2600); // (128-1)/2.
      DirectionCoordinate refCoord=dirCoords;
      refCoord.setIncrement(Vector<Double>(2, cell.get("rad").getValue()));
      refCoord.setReferencePixel(Vector<Double>(2,512.0));
      /*DirectionCoordinate dirCoordsBig(
          MDirection::J2000, Projection(Projection::SIN), 135 * C::pi / 180.0,
          60 * C::pi / 180.0, -5. * C::pi / 180.0 / 3600.0,
          5. * C::pi / 180.0 / 3600.0, xform, 63.5, 63.5);*/
      DirectionCoordinate dirCoordsSmall(
          MDirection::J2000, Projection(Projection::SIN), RA.get("rad").getValue(),
          DEC.get("rad").getValue(), Quantity(-36.0, "arcsec").get("rad").getValue(),
          Quantity(36.0, "arcsec").get("rad").getValue(), xform, 260, 260);
      //Vector<String> units(2);
      //units = "deg";
      //dirCoords.setWorldAxisUnits(units);

      // StokesCoordinate
      Vector<Int> iquv(1);
      iquv(0) = Stokes::I;
      StokesCoordinate stokesCoordsBad(iquv);
      Vector<Int> stoks(4);
      stoks(0) = Stokes::RR;
      stoks(1) = Stokes::RL;
      stoks(2) = Stokes::LR;
      stoks(3) = Stokes::LL;
      StokesCoordinate stokesCoordsGood(stoks);

      // SpectralCoordinate
      
      SpectralCoordinate spectralCoords3(MFrequency::TOPO, leFreq,
                                        128 * 1.0E+6, 0, 2000.40575 * 1.0E+6);
      //units.resize(1);
      //units = "Hz";
      spectralCoords3.setWorldAxisUnits(Vector<String>(1,"Hz"));

      SpectralCoordinate spectralCoords(MFrequency::TOPO, 1. * 1.0E+9,
                                         1 * 1.0E+9, 0, 1. * 1.0E+9);
      //spectralCoords.setWorldAxisUnits(units);
      SpectralCoordinate spectralCoords4(MFrequency::TOPO, 6.035328 * 1.0E+9,
                                         20 * 1.0E+3, 0, 6.035328 * 1.0E+9);

      
      //spectralCoords3.setWorldAxisUnits(units);

      // CoordinateSystem
      coordsys.addCoordinate(refCoord);
      coordsys.addCoordinate(stokesCoordsGood);
      coordsys.addCoordinate(spectralCoords);
      coordsys3.addCoordinate(dirCoords);
      coordsys3.addCoordinate(stokesCoordsGood);
      coordsys3.addCoordinate(spectralCoords3);
      coordsys3Small.addCoordinate(dirCoordsSmall);
      coordsys3Small.addCoordinate(stokesCoordsGood);
      coordsys3Small.addCoordinate(spectralCoords3);
      /*coordsys3Big.addCoordinate(dirCoordsBig);
      coordsys3Big.addCoordinate(stokesCoordsGood);
      coordsys3Big.addCoordinate(spectralCoords3);

      coordsys3Small.addCoordinate(dirCoordsSmall);
      coordsys3Small.addCoordinate(stokesCoordsGood);
      coordsys3Small.addCoordinate(spectralCoords3);
*/
      coordsys4.addCoordinate(dirCoords);
      coordsys4.addCoordinate(stokesCoordsGood);
      coordsys4.addCoordinate(spectralCoords4);
    }

    
    cerr << "Size " << findConvSize(coordsys3) << endl;
    cerr << "Size optim " << findConvSize(coordsys) << endl;
    /*String name("tab1");*/
    TiledShape ts(IPosition(4, 512, 512, 4, 1));
    TiledShape ts2(IPosition(4, 5200, 5200, 4, 1));
    TiledShape ts2Small(IPosition(4, 520, 520, 4, 1));
    //PagedImage<Complex> im1(ts, coordsys, name);
    //PagedImage<Float> im2(ts, coordsys, "tab2");
    PagedImage<Complex> im3(ts, coordsys, "tab3");
    //PagedImage<Float> im4(ts2, coordsys3, "tab4");

    //im1.set(Complex(1.0, 1.0));
    //im2.set(0.0);
    im3.set(Complex(1.0));
    //im4.set(0.0);
    PagedImage<Complex> diagIm(ts2, coordsys3, "DiagMuellerSquint.im");
    PagedImage<Complex> origIm(ts2, coordsys3, "OrigMuellerSquint.im");
    PagedImage<Complex> origNo(ts2, coordsys3, "OrigMuellerNoSquint.im");
    PagedImage<Complex> avgPB(ts2, coordsys3, "DiagMuellerNoSquint.im");
    PagedImage<Complex> diagImSmall(ts2Small, coordsys3Small, "DiagMuellerSquintSmall.im");
    diagIm.set(Complex(1.0));
    diagImSmall.set(1.0);
    diagImSmall.set(Complex(1.0));
    avgPB.set(1.0);
    ///////////////////////////////////////////

    { // begin tests

      
      ////////////////////////////////////////////////

      {
        // const char *sep=" ";
        // cerr << "getenv " << getenv("CASAPATH") << endl;
        // char *aipsPath = strtok(getenv("CASAPATH"),sep);

        // if (aipsPath == NULL)
        //   throw(AipsError("CASAPATH not found."));

        String msFileName;
        msFileName = "/home/casa/data/master/regression/unittest/concat/input/"
                     "A2256LC2_4.5s-1.ms";

        cout << "Reading " << msFileName << endl;

        MS ms(msFileName, Table::Old);

        Block<int> sort(4);
        sort[2] = MS::FIELD_ID;
        sort[3] = MS::ARRAY_ID;
        sort[1] = MS::DATA_DESC_ID;
        sort[0] = MS::TIME;
        vi::VisibilityIterator2 vi2(ms, vi::SortColumns(), false);
        vi::VisBuffer2 *vb = vi2.getVisBuffer();
        // ROVisibilityIterator vi(ms, sort, 0.);

        // VisBuffer vb(vi);

        cout << "new EVLAAperture" << endl;
        refim::EVLAAperture apB;
        apB.cacheVBInfo("EVLA", 25.0);
        Int count = 0;
        count = 0;
        vi2.originChunks();
        vi2.origin();
        //for (vi2.originChunks(); vi2.moreChunks(); vi2.nextChunk())
        {
          cout << "next chunk" << endl;
          //for (vi2.origin(); vi2.more(); vi2.next()) 
          {
            // the following can't be tested here since getVisParams is
            // protected
            // cout << "band id " << apB.getVisParams(vb) << endl;

            Int cfKey = 0;
          auto aptr=std::make_shared<refim::EVLAAperture>();
          std::shared_ptr<refim::WPConvFunc>wptr;
          aptr->cacheVBInfo("EVLA", 25.0);
          refim::AWConvFunc a(aptr, wptr);
          Array<Complex>  con;
          Array<Complex>  wtcon;
          Vector<Int> sup;
          Vector<Double> freqs={2e9, 2.2e9, 2.4e9};
          Timer timo;
          timo.mark();
          Int npix=-1;
          a.makeAConvFunc(con, wtcon, coordsys3, sup, npix, freqs, False, 0.0);
          cerr << "Shapes " << con.shape() << max(con) << " pixels used" << npix<< endl;
          timo.show("@@@makeAconv");
          {
            PagedImage<Complex> loo(con.shape(), coordsys3, "ELCONV");
            loo.put(con);
            PagedImage<Complex> loo2(con.shape(), coordsys3, "ELwtCONV");
            loo2.put(wtcon);
          }
           /* try {
              apB.applySky(im1, *vb, true, cfKey);
            } catch (AipsError x) {
              cout << "Caught expected error: " << x.getMesg() << endl;
            }*/
            cerr << "im3 shape " << im3.shape() << endl;
            Double pa=vb->feedPa(vb->time()(0))(0);
            pa=0.0;
            Timer tim;
            tim.mark();
            
            apB.applySky(im3, pa, true, cfKey, 0, leFreq);
            tim.show("applySky M0");
            tim.mark();
            FFT2D ftim3;
            ftim3.c2cFFT(im3);
            Float factor=cell.get("rad").getValue()/imcell.get("rad").getValue();
            //resampleCopy(origIm, im3, factor, factor);
            Float newfactor=1.0/(factor*512.0/5200);
            fakeConv(im3, origIm, newfactor, 30);
            Int newImsize= Int(std::round(512.0*newfactor*2.0))/2;
            cerr << "NEWIMSize " << newImsize << " newfactor "<< newfactor <<  endl;
            CoordinateSystem FTOrigcs=refim::SynthesisUtils::makeUVCoords(coordsys3, diagIm.shape());
            CoordinateSystem FTconvcs=refim::SynthesisUtils::makeUVCoords(coordsys, im3.shape());
            cerr << " orig ft incr " << FTOrigcs.increment() << "  " << FTconvcs.increment() << endl;
            tim.show("resample time");
            /*IPosition blc(4, 0, 0, 0, 0);
            IPosition trc=im3.shape()-1;
            trc[2]=0;
            Slicer sl(blc, trc, Slicer::endIsLast);
            SubImage<Complex> inplane(im3, sl, false);
            SubImage<Complex> outplane(origIm, sl, true);
            outplane.copyData(inplane);
            im3.set(1.0);
            tim.mark();
            apB.applySky(im3, pa, true, cfKey, 5, leFreq);
            tim.show("applySky M5");
            blc[2]=1;
            trc[2]=1;
            sl=Slicer(blc, trc, Slicer::endIsLast);
            outplane=SubImage<Complex>(origIm, sl, true);
            cerr << "outshape " << outplane.shape() << endl;
            outplane.copyData(inplane);
            im3.set(1.0);
            tim.mark();
            apB.applySky(im3, pa, true, cfKey, 10, leFreq);
            tim.show("applySky M10");
            blc[2]=trc[2]=2;
            sl=Slicer(blc, trc, Slicer::endIsLast);
            outplane=SubImage<Complex>(origIm, sl, true);
            outplane.copyData(inplane);
            im3.set(1.0);
            tim.mark();
            apB.applySky(im3, pa, true, cfKey, 15, leFreq);
            tim.show("applySky M15");
            blc[2]=trc[2]=3;
            sl=Slicer(blc, trc, Slicer::endIsLast);
            outplane=SubImage<Complex>(origIm, sl, true);
            outplane.copyData(inplane);
            */
            Double tbig=0.0;
            Double tsmall=0.0;
            tim.mark();
            for (int k=0; k<1; ++k){
              diagIm.set(1.0);
              pa=Double(k)*20.0*C::pi/180.0;
              aa.applyDiagSkyJones(diagIm, pa);
              PagedImage<Complex> outim(diagIm.shape(), diagIm.coordinates(), "Beam_"+String::toString(k*20)+".im");
              outim.copyData(diagIm);
            }
              tbig+=tim.real();
              tim.show("applyDiagSkyJones M0,5,10,15 ");
            tim.mark();
            aa.applyDiagSkyJones(diagImSmall, pa);
            tsmall += tim.real();
            tim.show("applyDiagSkyJonesSmall M0,5,10,15 ");
            CoordinateSystem FTbigcs=diagIm.coordinates();
            CoordinateSystem FTSmcs=diagImSmall.coordinates();
            FTbigcs=refim::SynthesisUtils::makeUVCoords(FTbigcs, diagIm.shape());
            FTSmcs=refim::SynthesisUtils::makeUVCoords(FTSmcs, diagImSmall.shape());
            PagedImage<Complex> bigConv(diagIm.shape(), FTbigcs, "BigConv.im");
            PagedImage<Complex> smConv(diagImSmall.shape(), FTSmcs, "SmallConv.im");
            tim.mark();
            bigConv.copyData(diagIm);
            tbig+=tim.real();
            tim.mark();
            smConv.copyData(diagImSmall);
            tsmall += tim.real();
            FFT2D ftbig;
            FFT2D ftsm;
            tim.mark();
            ftbig.c2cFFT(bigConv);
            tbig += tim.real();
            tim.show("FFT big M0,5,10,15 ");
            tim.mark();
            ftsm.c2cFFT(smConv);
            tsmall += tim.real();
            tim.show("FFT small M0,5,10,15 ");
            cerr << "Time taken big= " << tbig << " small= " << tsmall << endl; 
            
            ///No squint now
            /*m3.set(1.0);
            tim.mark();
            apB.applySky(im3, pa, false, cfKey, 0, leFreq);
            tim.show("applySkyNo M0");
            blc[2]=0;
            trc[2]=0;
            sl=Slicer(blc, trc, Slicer::endIsLast);
            outplane=SubImage<Complex>(origNo, sl, true);
            outplane.copyData(inplane);
            im3.set(1.0);
            tim.mark();
            apB.applySky(im3, pa, false, cfKey, 5, leFreq);
            tim.show("applySkyNo M5");
            blc[2]=1;
            trc[2]=1;
            sl=Slicer(blc, trc, Slicer::endIsLast);
            outplane=SubImage<Complex>(origNo, sl, true);
            cerr << "outshape " << outplane.shape() << endl;
            outplane.copyData(inplane);
            im3.set(1.0);
            tim.mark();
            apB.applySky(im3, pa, false, cfKey, 10, leFreq);
            tim.show("applySkyNo M10");
            blc[2]=trc[2]=2;
            sl=Slicer(blc, trc, Slicer::endIsLast);
            outplane=SubImage<Complex>(origNo, sl, true);
            outplane.copyData(inplane);
            im3.set(1.0);
            tim.mark();
            apB.applySky(im3, pa, false, cfKey, 15, leFreq);
            tim.show("applySkyNo M15");
            blc[2]=trc[2]=3;
            sl=Slicer(blc, trc, Slicer::endIsLast);
            outplane=SubImage<Complex>(origNo, sl, true);
            outplane.copyData(inplane);
            */
            
            tim.mark();
            aa.applyAvgSkyJones(avgPB);
            tim.show("applyAvgPB M0,5,10,15 ");
          }
          count++;
          cerr << "COUNT " << count << endl;
          // five rounds is enough
          //if (count > 5)
          //  break;
        }
        return 0;
        cout << endl
             << "******** ApplySky to unity image with squint and without "
                "***********"
             << endl
             << endl;

        count = 0;
        Int count2 = 0;

        for (vi2.originChunks(); vi2.moreChunks(); vi2.nextChunk()) {
          cout << "next chunk" << endl;
          for (vi2.origin(); vi2.more(); vi2.next()) {

            if (count2 == 0) { // first occurence

              Int cfKey = 0;

              PagedImage<Complex> im5(ts2, coordsys3, "pb_squintEVLA");
              im5.set(Complex(1.0, 1.0));
              apB.applySky(im5, *vb, true, cfKey);

              PagedImage<Complex> im6(ts2, coordsys3, "pb_nosquintEVLA");
              im6.set(Complex(1.0, 1.0));
              apB.applySky(im6, *vb, false, cfKey);

              // 	      PagedImage<Complex> im13(ts2, coordsys3Big,
              // "pb_squintEVLABig"); 	      im13.set(Complex(1.0,1.0));
              // 	      apB.applySky(im13, vb, true, cfKey);

              // 	      PagedImage<Complex> im14(ts2, coordsys3Small,
              // "pb_squintEVLASmall"); im14.set(Complex(1.0,1.0));
              // 	      apB.applySky(im14, vb, true, cfKey);

              count2++;

            } else if (count2 > 100) { // pick later occurance (to vary PA)

              Int cfKey = 0;

              PagedImage<Complex> im15(ts2, coordsys3, "pb2_squintEVLA");
              im15.set(Complex(1.0, 1.0));
              apB.applySky(im15, *vb, true, cfKey);

              // 	      PagedImage<Complex> im17(ts2, coordsys3Big,
              // "pb2_squintEVLABig"); 	      im17.set(Complex(1.0,1.0));
              // 	      apB.applySky(im17, vb, true, cfKey);

              // 	      PagedImage<Complex> im18(ts2, coordsys3Small,
              // "pb2_squintEVLASmall"); im18.set(Complex(1.0,1.0));
              // 	      apB.applySky(im18, vb, true, cfKey);

              // 	      PagedImage<Float> im18b(ts2, coordsys3Small,
              // "pb2_squintEVLASmall_float"); 	      im18b.set(1.0);
              // 	      apB.applySky(im18b, vb, true, cfKey);

              count2++;

              count = -1;
            } else {
              count2++;
            }
          }
          if (count < 0)
            break;
          count++;
        }
      }
      ////////////////////////////////////////////////

      {
        // const char *sep=" ";
        // char *aipsPath = strtok(getenv("CASAPATH"),sep);
        // if (aipsPath == NULL)
        //   throw(AipsError("CASAPATH not found."));

        String msFileName;
        msFileName = "/tmp/evla-highres-sample.ms";
        cout << "Reading " << msFileName << endl;

        MS ms(msFileName, Table::Old);
        vi::VisibilityIterator2 vi2(ms, vi::SortColumns(), false);
        vi::VisBuffer2 *vb = vi2.getVisBuffer();
        /*Block<int> sort(4);
        sort[2] = MS::FIELD_ID;
        sort[3] = MS::ARRAY_ID;
        sort[1] = MS::DATA_DESC_ID;
        sort[0] = MS::TIME;

        ROVisibilityIterator vi(ms, sort, 0.);

        VisBuffer vb(vi);*/

        cout << "new EVLAAperture" << endl;
        refim::EVLAAperture apB;

        cout << endl
             << "******** ApplySky to unity image with squint and without 2"
                "***********"
             << endl
             << endl;

        Int count = 0;
        Int count2 = 0;

        for (vi2.originChunks(); vi2.moreChunks(); vi2.nextChunk()) {
          cout << "next chunk" << endl;
          for (vi2.origin(); vi2.more(); vi2.next()) {

            if (count2 == 0) { // first occurence

              Int cfKey = 0;

              PagedImage<Complex> im5(ts2, coordsys4, "pb3_squintEVLA");
              im5.set(Complex(1.0, 1.0));
              apB.applySky(im5, *vb, true, cfKey);

              PagedImage<Complex> im6(ts2, coordsys4, "pb3_nosquintEVLA");
              im6.set(Complex(1.0, 1.0));
              apB.applySky(im6, *vb, false, cfKey);

              count2++;

            } else if (count2 > 100) { // pick later occurance (to vary PA)

              Int cfKey = 0;

              PagedImage<Complex> im15(ts2, coordsys4, "pb4_squintEVLA");
              im15.set(Complex(1.0, 1.0));
              apB.applySky(im15, *vb, true, cfKey);

              count2++;

              count = -1;
            } else {
              count2++;
            }
          }
          if (count < 0)
            break;
          count++;
        }
      }

    } // end tests

  } catch (AipsError x) {
    cout << "Caught Error: " << x.getMesg() << endl;
    exit(1);
  }

  cout << "OK" << endl;
  exit(0);
}
