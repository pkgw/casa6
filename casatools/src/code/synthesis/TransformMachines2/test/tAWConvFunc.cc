//# tAWConvFunc.cc: This program tests the AWConvFunc class
//# Copyright (C) 2023
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
#include <synthesis/TransformMachines2/AWConvFuncHolder.h>
#include <synthesis/TransformMachines2/Utils.h>
#include <synthesis/TransformMachines/PBMath1DAiry.h>
#include <synthesis/ImagerObjects/SynthesisUtilMethods.h>

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

int mainDoofus() {
  try {

    // Form coordinate systems and test images
    CoordinateSystem coordsys;
    CoordinateSystem coordsys3;
    CoordinateSystem coordWsys;
   
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
 
      //Vector<String> units(2);
      //units = "deg";
      //dirCoords.setWorldAxisUnits(units);

      // StokesCoordinate
     
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

      SpectralCoordinate spectralCoords(MFrequency::TOPO, leFreq,
                                         128 * 1.0E+6, 0,  2000.40575* 1.0E+6);
   
      //spectralCoords3.setWorldAxisUnits(units);

      // CoordinateSystem
      coordsys.addCoordinate(refCoord);
      coordsys.addCoordinate(stokesCoordsGood);
      coordsys.addCoordinate(spectralCoords3);
      coordsys3.addCoordinate(dirCoords);
      coordsys3.addCoordinate(stokesCoordsGood);
      coordsys3.addCoordinate(spectralCoords3);
      coordWsys.addCoordinate(refCoord);
      coordWsys.addCoordinate(spectralCoords);
     


    }

    
    cerr << "Size " << findConvSize(coordsys3) << endl;
    cerr << "Size optim " << findConvSize(coordsys) << endl;
    /*String name("tab1");*/
    TiledShape ts(IPosition(4, 512, 512, 4, 1));
    TiledShape ts2(IPosition(4, 5200, 5200, 4, 1));
    PagedImage<Complex> diagIm(ts, coordsys, "PB.im");
    PagedImage<Complex> im3(ts, coordsys, "Wterm.im");
    
    //PagedImage<Float> im4(ts2, coordsys3, "tab4");
    FFT2D ftim;
    //im1.set(Complex(1.0, 1.0));
    //im2.set(0.0);
    im3.set(Complex(1.0));
    //im4.set(0.0);
    //PagedImage<Complex> diagIm(ts2, coordsys3, "DiagMuellerSquint.im");
    //PagedImage<Complex> origIm(ts2, coordsys3, "OrigMuellerSquint.im");
    //PagedImage<Complex> origNo(ts2, coordsys3, "OrigMuellerNoSquint.im");
    //PagedImage<Complex> avgPB(ts2, coordsys3, "DiagMuellerNoSquint.im");
    
    //diagIm.set(Complex(1.0));
   
   // avgPB.set(1.0);
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
         Timer timo;
        Block<int> sort(4);
        sort[2] = MS::FIELD_ID;
        sort[3] = MS::ARRAY_ID;
        sort[1] = MS::DATA_DESC_ID;
        sort[0] = MS::TIME;
        vi::VisibilityIterator2 vi2(ms, vi::SortColumns(), false);
        vi::VisBuffer2 *vb = vi2.getVisBuffer();     

      
        vi2.originChunks();
        vi2.origin();
        //for (vi2.originChunks(); vi2.moreChunks(); vi2.nextChunk())
        
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
          diagIm.set(1.0);
          aptr->applyDiagSkyJones(diagIm,0.0);
          refim::AWConvFunc a(aptr, wptr);
          Array<Complex>  con; 
          Array<Complex>  wtcon;
          Vector<Int> sup;
          Vector<Double> freqs(32);
          for( int k=0; k < 32; ++k){
              freqs(k)=2.0000001e9 + Double(k)/32.0*2.0e9;
            
          }
         
          timo.mark();
          Int npix;
          a.makeAConvFunc(con, wtcon, coordsys3, sup, npix, freqs, False, 0.0);
          cerr << "Shapes " << con.shape() << max(con) << endl;
          timo.show("@@@makeAconv");
          {
            PagedImage<Complex> loo(con.shape(), coordsys3, "ELCONV");
            loo.put(con);
            PagedImage<Complex> loo2(con.shape(), coordsys3, "ELwtCONV");
            loo2.put(wtcon);
          }
          
          refim::WPConvFunc wpc;
          Matrix<Complex> arr(512,512);
          timo.mark();
          wpc.makeSkyWFunc(arr, coordsys, 512, Double(2e5));
          timo.show("@@@ one W ");
          im3.putSlice(arr, IPosition(4,0));
          Vector<Double> wVals={0, 1e4, 3e4, 2e5};
          Cube<Complex> wCon;
          Vector<int> sups;
          timo.mark();
          wpc.makeWConvFuncs(wCon, sups, coordsys, 512, wVals);
          cerr << "supports " << sups << endl;
          timo.show("@@@ wconvs ");
           {
            PagedImage<Complex> loo(wCon.shape(), coordWsys, "ELWConv");
            loo.put(wCon);
           }
           Array<Complex> aWConv;
           Array<Complex> aWwtconv;
           Matrix<Int> awSupport;
           timo.mark();
           for (int k=0; k<1; ++k){
              Double pa=Double(k)*20.0*C::pi/180.0;
              a.makeAWConvFunc(aWConv, aWwtconv,coordsys,awSupport, npix, freqs, wVals, True, pa);
           }
          timo.show("AWCONF ");
          Vector<Double>pixW(wVals.nelements());
          indgen(pixW);
           CoordinateSystem fiveAxis=coordsys;
          TabularCoordinate tab(pixW, wVals, "m", "W");
          fiveAxis.addCoordinate(tab);
          PagedImage<Complex> noo(aWConv.shape(), fiveAxis, "AWConvVals");
          noo.put(aWConv);
        
          
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

void makeImCoordsys(CoordinateSystem& cs,  const String& msname ) {
  
  
    Record msrec;
    msrec.define("msname", msname);
    msrec.define("field", "0");
    msrec.define("spw", "0:5~51");
    Record parsrec;
    parsrec.defineRecord("ms0", msrec);
    //Vector<Double> start(20);
    //Vector<Double> end(20);
    //start(0)=1.412787e9;
    //end(0)=1.412787e9+2.5e4;
    //for (Int k=1; k < 20; ++k){
    //  start(k)=start(k-1)+2.5e4;
    //  end(k)=end(k-1)+2.5e4;
    //}
	  
    SynthesisParamsImage impars;
    Vector<Int> ims(2);ims[0]=5000; ims[1]=5000;
    impars.imsize=ims;
    Vector<Quantity> cells(2); cells[0]=Quantity(6, "arcsec"), cells[1]=Quantity(6,"arcsec");
    impars.cellsize=cells;
    impars.stokes="IQUV";
    Quantity RA; Quantity::read(RA, "18:00:45");
      Quantity DEC; Quantity::read(DEC,"78.28.04");
    impars.phaseCenter=MDirection(RA, DEC, MDirection::J2000);
    impars.nchan=20;
    //impars.freqStart=freqStart;
    //impars.freqStep=freqStep;
    impars.restFreq=Quantity(1.420, "GHz");
    //impars.nTaylorTerms=nTaylorTerms;
    //impars.refFreq=refFreq;
    impars.projection=Projection::SIN;
    impars.freqFrame=MFrequency::LSRK;
    //impars.overwrite=overwrite;
    //impars.startModel=startmodel;
    impars.velStart=Quantity(1320.0, "km/s");
    impars.velStep=Quantity(0.3, "km/s");
    //impars.veltype="OPTICAL";
    impars.veltype="RADIO";

    MeasurementSet ms(msname, Table::Old);
  Vector<Int>spwid = {0, 1};
    cs=impars.buildCoordinateSystemCore(ms, spwid,  0,  1000.0,  1200.0,  1010.0, 1190.0 );

  
  
  
}

int main(){
  String msFileName;
        msFileName = "/home/casa/data/master/regression/unittest/concat/input/"
                     "A2256LC2_4.5s-1.ms";

        cout << "Reading " << msFileName << endl;

        MS ms(msFileName, Table::Old);
         Timer timo;
         std::vector<Double> usefulFreqs;
        vi::VisibilityIterator2 vi2(ms, vi::SortColumns(), false);
        vi::VisBuffer2 *vb = vi2.getVisBuffer();     
        for (vi2.originChunks(); vi2.moreChunks(); vi2.nextChunk()) {
          for (vi2.origin(); vi2.more(); vi2.next()) {
              for (int somechan = 0; somechan < 5; ++somechan) {
                 usefulFreqs.push_back(vb->getFrequency(0,  somechan*5));
              }
            
          }
        }
        vi2.originChunks();
        vi2.origin();
        String observatory=(vb->subtableColumns().observation()).telescopeName()(0);
        cerr <<  "FReqs" <<  usefulFreqs <<  endl;
        std::sort(usefulFreqs.begin(),  usefulFreqs.end());
        auto last = std::unique(usefulFreqs.begin(),  usefulFreqs.end());
        usefulFreqs.erase(last,  usefulFreqs.end());
        cerr <<  "Freqs " <<  usefulFreqs <<  endl;
        Vector<Double>freqs(usefulFreqs);
        Vector<Double> wVals = {0, 1e3, 2e3,  4e3,  8e3,  1.2e4, 3e4};
        CoordinateSystem csys;
        makeImCoordsys(csys, msFileName);
       refim::AWConvFuncHolder awh(csys, 4000, 4000, False, 2.0*C::pi/18.0,  observatory,  4);
       timo.mark();
       awh.addConvFunc(freqs,  wVals,  0.0);
       /////Get and save convfunc
       {
        Array<Complex> convFunc(awh.getConvFunc());
        Vector<Double>pixW(convFunc.shape()[4]);
        indgen(pixW);
        CoordinateSystem fiveAxis=csys;
        TabularCoordinate tab(pixW, pixW, "m", "W");
        fiveAxis.addCoordinate(tab);
        PagedImage<Complex> noo(convFunc.shape(), fiveAxis, "AWConvVals");
        noo.put(convFunc);
       }
       /////////
       timo.show("@@Time taken for 8w *"+String::toString(freqs.nelements())+"freqs");
       Vector<Int> convPolMap; 
       Vector<Int> convChanMap;
       Vector<Int> convRowMap;
        for (vi2.originChunks(); vi2.moreChunks(); vi2.nextChunk()) {
          for (vi2.origin(); vi2.more(); vi2.next()) {
              awh.getConvIndices(convPolMap, convChanMap, convRowMap,*vb);
              cerr << "pol " << convPolMap << endl;
              cerr << "chan " << convChanMap << endl;
              cerr << "rowMap " << convRowMap << endl;
            
          }
        }
        
  
 return 0; 
}
