/*
 * tSynthesisImager.cc
 *SynthesisImager.cc: test of SynthesisImager
//# Copyright (C) 2021
//# Associated Universities, Inc. Washington DC, USA.
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
 *  Created on: April 01 2021
 *      
 */


#include <casa/iostream.h>
#include <casa/aips.h>
#include <casa/Exceptions/Error.h>
#include <casa/BasicSL/String.h>
#include <casa/Containers/Block.h>
#include <measures/Measures/MRadialVelocity.h>
#include <coordinates/Coordinates/CoordinateSystem.h>
#include <casa/Logging/LogIO.h>
#include <synthesis/ImagerObjects/SynthesisImagerVi2.h>
#include <synthesis/ImagerObjects/SynthesisUtilMethods.h>
#include <synthesis/ImagerObjects/SIImageStore.h>
#include <synthesis/ImagerObjects/SimpleSIImageStore.h>
#include <imageanalysis/Utilities/SpectralImageUtil.h>
#include <lattices/Lattices/LatticeConcat.h>
#include <images/Images/PagedImage.h>
#include <images/Images/ImageConcat.h>
#include <images/Images/SubImage.h>
#include <casa/namespace.h>
#include <images/Images/TempImage.h>
#include <coordinates/Coordinates/CoordinateUtil.h>
#include <ms/MSSel/MSSourceIndex.h>
#include <synthesis/TransformMachines2/test/MakeMS.h>

int main(int argc, char **argv)
{
  using namespace std;
using namespace casacore;
  using namespace casa;
using namespace casacore;
  using namespace casa::test;
  try{


    
    if (argc<2) {
      cout <<"Usage: tAWProjectCube ms-table-name   "<<endl;
      exit(1);
    }
    String msname=String(argv[1]);
    MeasurementSet thems(msname, TableLock::UserNoReadLocking, Table::Old);
    	 
	  
    /*MDirection thedir(Quantity(20.0, "deg"), Quantity(20.0, "deg"));
      String msname("Test2.ms");
      const Int numchan=10;
      MakeMS::makems(msname, thedir, 1.0e9, 1e8, numchan, 20);
      //MakeMS::makems(msname, thedir, 1.5e9, 1e6, numchan, 20);
      
      MeasurementSet thems(msname, TableLock::UserNoReadLocking, Table::Update);
      thems.markForDelete();
      thems.lock(True);
      MSColumns(thems).data().fillColumn(Matrix<Complex>(4,numchan, Complex(6.66e-2)));
      MSColumns(thems).correctedData().fillColumn(Matrix<Complex>(4,numchan, Complex(6.66e-2)));
      thems.flush();
    */  
   
    CountedPtr<SynthesisImager> imgr = new SynthesisImagerVi2();
    SynthesisParamsSelect selpars;
    selpars.msname=msname; selpars.spw="*"; selpars.freqframe=MFrequency::LSRK;
    selpars.field="*"; selpars.usescratch=false; selpars.readonly=true;
    imgr->selectData(selpars);
    cout <<"--Imager created for MeasurementSet object. " << endl;
    thems.unlock();
    MeasurementSet tab(msname);
    MDirection phasecenter=MSFieldColumns(tab.field()).phaseDirMeas(0,0.0);
    Quantity freqBeg=MSSpWindowColumns(tab.spectralWindow()).chanFreqQuant()(0)(IPosition(1,0));
    Int ndataChan=MSSpWindowColumns(tab.spectralWindow()).numChan()(0);
    Quantity freqWidth=MSSpWindowColumns(tab.spectralWindow()).chanFreqQuant()(0)(IPosition(1,ndataChan-1));
    freqWidth-=freqBeg;
    int nx,ny;
    nx = ny = 1024;
    Quantity cellx( 10.0, "arcsec" );
    Quantity celly( 10.0, "arcsec" );
    //Quantity cellx( 0.5, "arcsec" );
    //Quantity celly( 0.5, "arcsec" );
    //Vector<Int> spwids(2);
    String stokes="I";
    //Int nchan=1;
    cerr << "nx=" << nx << " ny=" << ny
         << " cellx='" << cellx.getValue() << cellx.getUnit()
         << "' celly='" << celly.getValue() << celly.getUnit()
         << "' spwids= *" 
         << " field=" <<   0 << endl;
	  ////lets do a cube of ndatachan
	  //nchan=ndataChan;
          
    freqWidth /= Double(ndataChan);
    String imageName("test_cube_image1");
    String cfCache("test_cube.cf");
    SynthesisParamsImage impars;
    SynthesisParamsGrid gridpars;
    Vector<Quantity> qCellSize(2); qCellSize[0]=cellx; qCellSize[1]=celly;
    Vector<Int> ims(2);ims[0]=nx, ims[1]=ny;
            // Set up the parameters common for all modes
    impars.cellsize=qCellSize;
    impars.imsize=ims;
    impars.phaseCenter=phasecenter;
    impars.stokes=stokes;
    impars.freqStart=freqBeg;
    impars.freqStep=freqWidth;
    impars.restFreq=Vector<Quantity>(1,Quantity(1.420, "GHz"));
    impars.nchan=200;
    impars.mode="cube";
    impars.imageName=imageName;
    gridpars.facets=1;
    //	  gridpars.imageName=imageName;
    gridpars.imageName=imageName;
    gridpars.aTermOn=true;
    gridpars.psTermOn=false;
    gridpars.mTermOn=true;
    gridpars.wbAWP=true;
    gridpars.gridder=String("awprojectft");
    gridpars.ftmachine=String("awprojectft");
    gridpars.wprojplanes=4;
    gridpars.cfCache=cfCache;
    gridpars.conjBeams=true;
    ////lets do a cube of ndatachan
    //nchan=ndataChan;
    //freqWidth /= Double(nchan);
            
    //impars.nchan=nchan;
    //impars.freqStep=freqWidth;
    Record normpars;
    normpars.define("pblimit", 0.1); normpars.define("nterms",1); normpars.define("facets",1);
    normpars.define("normtype", "flatnoise"); normpars.define("workdir",".");
    normpars.define("deconvolver","hogbom");normpars.define("imagename",imageName); normpars.define("restoringbeam","");
    normpars.define("psfcutoff",0.35);
    imgr->defineImage(impars,gridpars);
    imgr->normalizerinfo(normpars);
    imgr->weight("natural");
    Record rec;
    imgr->executeMajorCycle(rec);
    imgr->makePSF();
            //imgr->makePSF(useViVb2);
    CountedPtr<SIImageStore> images=imgr->imageStore(0);
    if(images.null())
      throw(AipsError("Did not get shared_ptr "));
    images->dividePSFByWeight();
    images->divideResidualByWeight();

    Int nchan=(imgr->imageStore(0))->getShape()(3);
    
    
    cerr << "---------------Done with one monolithic cube ========================" << endl;
    ////done with cube
    
    cerr << "==============starting imaging by slices ============================" << endl;
    //Cube slicing 
    {
      
      imgr = new SynthesisImagerVi2();
      impars.imageName=String("test_cubesliced_image2");
      //	    imgr->defineImage(impars,gridpars);
      imgr->selectData(selpars);
      gridpars.cfCache="test_cube_slice.cf";
      imgr->defineImage(impars, gridpars);
      CountedPtr<SIImageStore> si=imgr->imageStore(0);
      CountedPtr<ImageInterface<Float> > resid=si->residual();
      CountedPtr<ImageInterface<Float> > psf=si->psf();
      CountedPtr<ImageInterface<Float> > wgt=si->weight();
      CountedPtr<ImageInterface<Float> > sumwt=si->sumwt();
      CountedPtr<ImageInterface<Float> > mod=si->model();
      CountedPtr<ImageInterface<Float> > restor=si->image();
      CountedPtr<ImageInterface<Float> > pb=si->pb();
      /////Imaging a channel at a time
      ////you could do a chunk at a time if memory allows
      
      {
        Int chunksize=4;
        for (Int k=0; k < nchan/chunksize; ++k){
          Int startchan=k*chunksize;
          Int endchan=k*chunksize+chunksize-1;
          
          std::shared_ptr<ImageInterface<Float> >subresid(SpectralImageUtil::getChannel(*resid, startchan, endchan, true));
          std::shared_ptr<ImageInterface<Float> >subpsf(SpectralImageUtil::getChannel(*psf, startchan, endchan, true));
          std::shared_ptr<ImageInterface<Float> > subwgt(SpectralImageUtil::getChannel(*wgt, startchan, endchan, true));
          std::shared_ptr<ImageInterface<Float> > subsumwt(SpectralImageUtil::getChannel(*sumwt, startchan, endchan, true));
          std::shared_ptr<ImageInterface<Float> > submod( SpectralImageUtil::getChannel(*mod, startchan, endchan, true));
          std::shared_ptr<ImageInterface<Float> > subrestor(SpectralImageUtil::getChannel(*restor, startchan, endchan, true));
          std::shared_ptr<ImageInterface<Float> > subpb(SpectralImageUtil::getChannel(*pb, startchan, endchan, true));
          /*
            std::shared_ptr<ImageInterface<Float> >subresid=std::make_shared<SubImage<Float> >(SpectralImageUtil::getChannel(*resid, k, k, true));
            
            std::shared_ptr<ImageInterface<Float> >subpsf= std::make_shared<SubImage<Float> >(SpectralImageUtil::getChannel(*psf, k, k, true));
            std::shared_ptr<ImageInterface<Float> > subwgt= std::make_shared<SubImage<Float> >(SpectralImageUtil::getChannel(*wgt, k, k, true));
            std::shared_ptr<ImageInterface<Float> > subsumwt= std::make_shared<SubImage<Float> >(SpectralImageUtil::getChannel(*sumwt, k, k, true));
            std::shared_ptr<ImageInterface<Float> > submod=std::make_shared<SubImage<Float> >( SpectralImageUtil::getChannel(*mod, k, k, true));
            std::shared_ptr<ImageInterface<Float> > subrestor= std::make_shared<SubImage<Float> >(SpectralImageUtil::getChannel(*restor, k, k, true));
          */
          CountedPtr<SIImageStore> subImStor=new SimpleSIImageStore(submod, subresid, subpsf,  subwgt, subrestor, nullptr, subsumwt, nullptr, subpb, nullptr, nullptr, true);

          
          //String freqBeg=String::toString(SpectralImageUtil::worldFreq(subresid->coordinates(), Double(-0.5)))+"Hz";
          //String freqEnd=String::toString(SpectralImageUtil::worldFreq(subresid->coordinates(), Double(0.5)))+"Hz";
          //	    CountedPtr<SIImageStore> subImStor=new SIImageStore(submod, subresid, subpsf, subwgt, subrestor, nullptr, nullptr, resid->coordinates(), "");
          
          SynthesisImagerVi2 subImgr;
          subImgr.selectData(selpars);
          subImgr.defineImage(subImStor, impars, gridpars);
          subImgr.weight("natural");
          Record rec;
          subImgr.normalizerinfo(normpars);
          subImgr.tuneSelectData();
          subImgr.setCubeGridding(False);
          subImgr.executeMajorCycle(rec);
          subImgr.makePSF();
          subImStor->dividePSFByWeight();
          subImStor->divideResidualByWeight();
        }
        
      }


     

    
      
            //We can do the division at the end
      //    si->dividePSFByWeight();
      //si->divideResidualByWeight();
    LatticeExprNode LEN = max( *(si->residual()) );
    cerr << "Sliced cube Max of whole residual image " << LEN.getFloat() << endl;
    LatticeExprNode psfmax = max( *(si->psf()) );
    LatticeExprNode psfmin = min( *(si->psf()) );
    cerr <<"Sliced cube Min max of whole psf "<< psfmin.getFloat() << " " << psfmax.getFloat() << endl;

    }
            
    {
        
        //After Normalization
        LatticeExprNode LENMaxRes = max( *(images->residual()) );
        cerr << "Monolithic cube Max of residual=" << LENMaxRes.getFloat()  << endl;
        LatticeExprNode psfmax = max( *(images->psf()) );
        LatticeExprNode psfmin = min( *(images->psf()) );
        cerr <<"Monolithic cube Min max of psf "<< psfmin.getFloat() << " " << psfmax.getFloat() << endl;
      }
           

         
	 
 


  }catch( AipsError e ){
    cout << "Exception ocurred." << endl;
    cout << e.getMesg() << endl;
    exit(-1);
  }
  cout << "OK" << endl;
  return 0;
};
