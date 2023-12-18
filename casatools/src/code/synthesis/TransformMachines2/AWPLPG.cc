//# AWPLPG.h: Implementation for a CPU based gridder for A,W (LPG= LowPerformanceGridder)
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

#include <synthesis/TransformMachines2/AWPLPG.h>
#include <synthesis/TransformMachines2/AWConvFuncHolder.h>
#include <synthesis/TransformMachines2/HetArrayConvFunc.h>
#include <synthesis/TransformMachines2/SimplePBConvFunc.h>
#include <msvis/MSVis/VisBuffer2.h>
#include <msvis/MSVis/VisibilityIteratorImpl2.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Vector.h>



namespace casa { //# NAMESPACE CASA - BEGIN
namespace refim {//# namespace for wstatimaging refactor
using namespace casacore;
using namespace casa;
using namespace casa::refim;

  
  AWPLPG::AWPLPG(SkyJones* sj, const Int nw,  Bool dosquint, const Double painc, MPosition mloc, String stokes,  const Bool usezero, const Bool useDoublePrec,  const casacore::Bool usePointing): MosaicFTNew(sj, mloc, stokes, Long(1000000), 16, usezero, True, False, usePointing), doSquint_p(dosquint), paInc_p(painc), nw_p(nw) {

    useDoubleGrid_p=True; //We'll always use double prec for this grid

  }
  AWPLPG::AWPLPG(const AWPLPG& other) : MosaicFTNew(other)
  {
    operator=(other);
  }

 AWPLPG& AWPLPG::operator=(const AWPLPG& other) {
  if(this!=&other) {

    //Do the base parameters
    MosaicFTNew::operator=(other);
    awConvs_p = other.awConvs_p;
    doSquint_p = other.doSquint_p;
    paInc_p = other.paInc_p;
    nw_p = other.nw_p;
    }
    return *this;
    
  }
  refim::FTMachine* AWPLPG::cloneFTM(){
      return new AWPLPG(*this);
  }
void AWPLPG::init(const vi::VisBuffer2& vb){
 MosaicFTNew::init(vb);
  
 
  //oversample if image is small
  //But not more than 5000 pixels
 convSampling=(max(nx, ny) < 50) ? 100: Int(ceil(5000.0/max(nx, ny)));
  if(convSampling <4) 
    convSampling=4;
 // TESTOO
  //convSampling = 1;
  // TESTOO
  
  CoordinateSystem cs=image->coordinates();
  
    SpectralCoordinate spCS = cs.spectralCoordinate(cs.findCoordinate(Coordinate::SPECTRAL));
    double f1, f2;
    nchan = image->shape()(3);
    spCS.toWorld(f1, double(-0.5));
    spCS.toWorld(f2, double(nchan)-0.5);
    auto frange=std::make_pair(f1, f2);
    
  
  if(awConvs_p.use_count()==0){
     String observatory=(vb.subtableColumns().observation()).telescopeName()(0);
    awConvs_p=std::make_shared<AWConvFuncHolder>((*image).coordinates(), nx, ny, 
                   doSquint_p, paInc_p, observatory, convSampling);
    vi::VisibilityIterator2 *vi= const_cast<VisibilityIterator2 *>(vb.getVi());
    
    if(pbConvFunc_p.null())
      pbConvFunc_p=new HetArrayConvFunc();
    if(sj_p)
      pbConvFunc_p->setSkyJones(sj_p.get());
      
    if((pbConvFunc_p->getVBUtil()).null()){
      if(vbutil_p.null()){
        vbutil_p=new VisBufferUtil(vb);
      }
      pbConvFunc_p->setVBUtil(vbutil_p);
    }
    std::vector<Double> freqs;
    std::vector<Double> pAs;
    Double maxW=0.0;
    for (vi->originChunks(); vi->moreChunks(); vi->nextChunk()) {
          for (vi->origin(); vi->more(); vi->next()) {
              std::vector<Double> chunkfreq;
              pbConvFunc_p->findUsefulChannels(chunkfreq, vb, frange);
              //cerr <<  "chunkfreq " <<  chunkfreq <<  endl;
              if (chunkfreq.size() > 0) {
                //cerr << "SPW " << vb.spectralWindows()(0) << " freqs " << Vector<Double>(chunkfreq) << endl;
                std::move(chunkfreq.begin(), chunkfreq.end(), std::back_inserter(freqs));
                double maxfreqused = *(std::max_element(chunkfreq.begin(), chunkfreq.end()));
                if (nw_p > 1) {
                  // 	maxW=max(maxW, max(abs(vb.uvw().row(2)*max(vb.getFrequencies(0))))/C::c);
                  maxW = max(maxW, max(abs(vb.uvw().row(2) * maxfreqused)) / C::c);
                }
              }
              if (doSquint_p)
                pAs.push_back(getPA(vb));
              //if(nw_p > 1)
              //  	maxW=max(maxW, max(abs(vb.uvw().row(2)*max(vb.getFrequencies(0))))/C::c);
          }
    }
    
    //return vi to origin
    vi->originChunks(); vi->origin();
    
    std::sort(freqs.begin(),  freqs.end());
    auto last = std::unique(freqs.begin(),  freqs.end());
    freqs.erase(last,  freqs.end());
    
    Double paMax=0.0;
    if(pAs.size()>1){
      std::sort(pAs.begin(), pAs.end());
      last=std::unique(pAs.begin(), pAs.end());
      pAs.erase(last, pAs.end());
      if(pAs.size()==1)
	paMax=pAs[0];
      if(pAs.size() >1){
	auto [minpa, maxpa]=std::minmax_element(begin(pAs), end(pAs));
        paMax=max(abs(*minpa), *maxpa);
      }
    }
    
    cerr <<  "PAMax in data " <<  paMax <<  endl;
    if (nw_p == 0)
      nw_p = 1;
    Vector<Double> wVals(nw_p,0);
    if(nw_p >1){
      Double st=maxW/(Double(nw_p-1)*Double(nw_p-1));
      for (int k=0; k <nw_p; ++k)
        wVals[k]=Double(k*k)*st;
    }
    (*awConvs_p).addConvFunc(Vector<Double>(freqs), wVals, paMax);
    
  }
  
}
  
 void AWPLPG::findConvFunction(const ImageInterface<Complex>& iimage, const vi::VisBuffer2& vb, const Matrix<Double>& rotuvw ){
  //
  // pbConvFunc_p.phasegradient
    convFunc.resize();
    convFunc.assign(awConvs_p->getConvFunc());
 
    weightConvFunc_p.resize();
    weightConvFunc_p.assign(awConvs_p->getWeightConvFunc());
    /*{ 
      ////TESTOO
      IPosition elshp = convFunc.shape().getFirst(4);
      IPosition elblc(5, 0);
      
      IPosition eltrc = convFunc.shape()-1;
      elblc[4] = eltrc[4];
      CoordinateSystem csysA = iimage.coordinates();
      Vector<Int> stoks(4);
      stoks(0) = Stokes::RR;
      stoks(1) = Stokes::RL;
      stoks(2) = Stokes::LR;
      stoks(3) = Stokes::LL;
      StokesCoordinate stokey(stoks);
      csysA.replaceCoordinate(stokey, 1);
      PagedImage<Complex> lastplane(elshp,  csysA,  "COOBOO");
      lastplane.put(convFunc(elblc,  eltrc).nonDegenerate());
      PagedImage<Complex> lastplaneW(elshp,  csysA,  "WOOBOO");
      lastplaneW.put(weightConvFunc_p(elblc,  eltrc).nonDegenerate());
    //////
    } */  
    convSizePlanes_p.resize();
    convSizePlanes_p = awConvs_p->getConvSizes();
    convSupportPlanes_p.resize();
    convSupportPlanes_p = awConvs_p->getConvSupports();
    awConvs_p->getConvIndices(convPolMap_p,  convChanMap_p,  convRowMap_p,  vb, rotuvw);
    //cerr <<  "min max convrowmap " <<  min(convRowMap_p) <<  "  " <<  max(convRowMap_p) <<  " supp " <<   max(convSupportPlanes_p) <<  " csize " << max(convSizePlanes_p) <<  " convchanmap "<< min(convChanMap_p) <<  "    " << max(convChanMap_p) << " convsamp " << convSampling << endl;
    std::vector<Int> pmapused=convPolMap_p.tovector();
    {
      std::sort(pmapused.begin(),  pmapused.end());
      auto last = std::unique(pmapused.begin(),  pmapused.end());
      pmapused.erase(last,  pmapused.end());
    }
     std::vector<Int> cmapused=convChanMap_p.tovector();
    {
      std::sort(cmapused.begin(),  cmapused.end());
      auto last = std::unique(cmapused.begin(),  cmapused.end());
      cmapused.erase(last,  cmapused.end());
    }
     std::vector<Int> rmapused=abs(convRowMap_p).tovector();
    {
      std::sort(rmapused.begin(),  rmapused.end());
      auto last = std::unique(rmapused.begin(),  rmapused.end());
      rmapused.erase(last,  rmapused.end());
    }
    //cerr << "pmap " << Vector<Int>(pmapused) << " cmp " << Vector<Int>(cmapused) << " rmap " << Vector<Int>(rmapused) << endl;
    pbConvFunc_p->rephaseConvFunc(iimage, vb, convSampling,  convFunc, weightConvFunc_p, pmapused, cmapused, rmapused,  MVDirection(-(movingDirShift_p.getAngle())), fixMovingSource_p);
    convSupport =max(convSupportPlanes_p);
    convSize = max(convSizePlanes_p);
   
    
 }
 
  /////==============================================
  //// some fortran defn
#define NEED_UNDERSCORES
#if defined(NEED_UNDERSCORES)
#define sectgmosd3 sectgmosd3_
#define sectdmos3 sectdmos3_
#define gmoswgtd2 gmoswgtd2_
#define locuvw locuvw_
#endif

extern "C" { 
  void locuvw(const Double*, const Double*, const Double*, const Int*, const Double*, const Double*, const Int*, 
	      Int*, Int*, Complex*, const Int*, const Int*, const Double*);
  void gmoswgtd2(const Int*/*nvispol*/, const Int*/*nvischan*/,
		const Int*/*flag*/, const Int*/*rflag*/, const Float*/*weight*/, const Int*/*nrow*/, 
		const Int*/*nx*/, const Int*/*ny*/, const Int*/*npol*/, const Int*/*nchan*/, 
		const Int*/*support*/, const Int*/*convsize*/, const Int*/*sampling*/, 
		const Int*/*chanmap*/, const Int*/*polmap*/,
		DComplex* /*weightgrid*/, Double* /*sumwt*/, const Complex*/*convweight*/, const Int*/*convplanemap*/, 
		const Int*/*convchanmap*/,  const Int*/*convpolmap*/, 
		const Int*/*nconvplane*/, const Int*/*nconvchan*/, const Int*/*nconvpol*/, const Int*/*rbeg*/, 
		const Int*/*rend*/, const Int*/*loc*/, const Int*/*off*/, const Complex*/*phasor*/);


  void sectgmosd3(const Complex* /*values*/,
		  Int* /*nvispol*/, Int* /*nvischan*/,
		  Int* /*dopsf*/, const Int* /*flag*/, const Int* /*rflag*/, const Float* /*weight*/,
		  Int* /* nrow*/, DComplex* /*grid*/, Int* /*nx*/, Int* /*ny*/, Int * /*npol*/, Int * /*nchan  */,
		  const Int*/*support*/, Int*/*convsize*/, Int*/*sampling*/, const Complex*/*convfunc*/,
		  const Int*/*chanmap*/, const Int*/*polmap*/,
		  Double*/*sumwgt*/, const Int*/*convplanemap*/,
		  const Int*/*convchanmap*/, const Int*/*convpolmap*/, 
		  Int*/*nconvplane*/, Int*/*nconvchan*/, Int* /*nconvpol*/,
		  const Int*/*x0*/,const Int*/*y0*/, const Int*/*nxsub*/, const Int*/*nysub*/, const Int*/*rbeg*/, 
		  const Int* /*rend*/, const Int*/*loc*/, const Int* /*off*/, const Complex*/*phasor*/);     

  void sectdmos3(Complex*,
  	      Int*,
  	      Int*,
  	      const Int*,
  	      const Int*,
		 Int*,
  	      const Complex*,
  	      Int*,
  	      Int*,
  	      Int *,
		 Int *,
		 const Int*,    //support
  	      Int*,
  	      Int*,
  	      const Complex*,
  	      const Int*,
  	      const Int*,
  	      const Int*,
	      const  Int*, 
	      const Int*, 
	      Int*, Int*, Int*,
		 //rbeg
		 const Int*,
		 const Int*,
		 const Int*,
		 const Int*,
		 const Complex*);

	     

}











  //===================================================
/*  void AWPLPG::put(const vi::VisBuffer2& vb, Int row, Bool dopsf,
		   FTMachine::Type type)
{


  
  
  Timer tim;
  tim.mark();
 
  matchChannel(vb);
 

  //cerr << "CHANMAP " << chanMap << endl;
  //No point in reading data if its not matching in frequency
  if(max(chanMap)==-1)
    return;

  //const Matrix<Float> *imagingweight;
  //imagingweight=&(vb.imagingWeight());
  Matrix<Float> imagingweight;
  getImagingWeight(imagingweight, vb);

  if(dopsf) type=FTMachine::PSF;

  Cube<Complex> data;
  //Fortran gridder need the flag as ints 
  Cube<Int> flags;
  Matrix<Float> elWeight;
  interpolateFrequencyTogrid(vb, imagingweight,data, flags, elWeight, type);
  
 

  Bool iswgtCopy;
  const Float *wgtStorage;
  wgtStorage=elWeight.getStorage(iswgtCopy);


  

  Bool isCopy;
  const Complex *datStorage=0;

  // cerr << "dopsf " << dopsf << " isWeightCopy " << iswgtCopy << "  " << wgtStorage<< endl;
  if(!dopsf)
    datStorage=data.getStorage(isCopy);
    
  
  // If row is -1 then we pass through all rows
  Int startRow, endRow, nRow;
  if (row==-1) {
    nRow=vb.nRows();
    startRow=0;
    endRow=nRow-1;
  } else {
    nRow=1;
    startRow=row;
    endRow=row;
  }
  
  // Get the uvws in a form that Fortran can use and do that
  // necessary phase rotation. 
  Matrix<Double> uvw(negateUV(vb));
  Vector<Double> dphase(vb.nRows());
  dphase=0.0;
 
  doUVWRotation_p=true;
  girarUVW(uvw, dphase, vb);
  refocus(uvw, vb.antenna1(), vb.antenna2(), dphase, vb);
  // This needs to be after the interp to get the interpolated channels
  //Also has to be after rotateuvw in case tracking is on
  findConvFunction(*image, vb);
  //cerr << "Put convsup " << convSupport << " max min convFunc " << max(convFunc) << "   " << min(convFunc) << "  "  << max(weightConvFunc_p) << min(weightConvFunc_p)  << "SHP " << convFunc.shape() << "   " << weightConvFunc_p.shape() << endl;
  //cerr << "convRowMap " << convRowMap_p  << " " << convChanMap_p << "  " << convPolMap_p << endl; 
  //nothing to grid here as the pointing resulted in a zero support convfunc
  if(convSupport <= 0)
    return;
  
  // Get the pointing positions. This can easily consume a lot 
  // of time thus we are for now assuming a field per 
  // vb chunk...need to change that accordingly if we start using
  // multiple pointings per vb.
  //Warning 

  // Take care of translation of Bools to Integer
  Int idopsf=0;
  if(dopsf) idopsf=1;
  
  
  Vector<Int> rowFlags(vb.nRows());
  rowFlags=0;
  rowFlags(vb.flagRow())=true;
  if(!usezero_p) {
    for (Int rownr=startRow; rownr<=endRow; rownr++) {
      if(vb.antenna1()(rownr)==vb.antenna2()(rownr)) rowFlags(rownr)=1;
    }
  }
  
  

  //cerr << "convSamp " << convSampling << " convsupp " << convSupport << " consize " << convSize << " convFunc " << convFunc.shape() << endl;
  //TESTOO

  //TESTOO
  
  //Tell the gridder to grid the weights too ...need to do that once only
  //Int doWeightGridding=1;
  //if(doneWeightImage_p)
  //  doWeightGridding=-1;
  Bool del;
  //    IPosition s(flags.shape());
  const IPosition& fs=flags.shape();
  //cerr << "flags shape " << fs << endl;
  std::vector<Int>s(fs.begin(), fs.end());
  Int nvp=s[0];
  Int nvc=s[1];
  Int nvisrow=s[2];
  Int csamp=convSampling;
  Bool uvwcopy; 
  const Double *uvwstor=uvw.getStorage(uvwcopy);
  Bool gridcopy;
  Bool convcopy;
  Bool wconvcopy;
  const Complex *convstor=convFunc.getStorage(convcopy);
  const Complex *wconvstor=weightConvFunc_p.getStorage(wconvcopy);
  Int nPolConv=convFunc.shape()[2];
  Int nChanConv=convFunc.shape()[3];
  Int nConvFunc=convFunc.shape()(4);
  Bool weightcopy;
  ////////
  Cube<Int> loc(2, nvc, nRow);
  Cube<Int> off(2, nvc, nRow);
  Matrix<Complex> phasor(nvc, nRow);
  Bool delphase;
  Complex * phasorstor=phasor.getStorage(delphase);
  const Double * visfreqstor=interpVisFreq_p.getStorage(del);
  const Double * scalestor=uvScale.getStorage(del);
  const Double * offsetstor=uvOffset.getStorage(del);
  Int * locstor=loc.getStorage(del);
  Int * offstor=off.getStorage(del);
  const Double *dpstor=dphase.getStorage(del);
  Int irow;
  Int nth=1;
#ifdef _OPENMP
  if(numthreads_p >0){
    nth=min(numthreads_p, omp_get_max_threads());
  }
  else{   
    nth= omp_get_max_threads();
  }
  //nth=min(4,nth);
#endif
  Double cinv=Double(1.0)/C::c;
 
  Int dow=0;
#pragma omp parallel default(none) private(irow) firstprivate(visfreqstor, nvc, scalestor, offsetstor, csamp, phasorstor, uvwstor, locstor, offstor, dpstor, dow, cinv) shared(startRow, endRow) num_threads(nth)  
{
#pragma omp for
  for (irow=startRow; irow<=endRow;irow++){
 
    locuvw(uvwstor, dpstor, visfreqstor, &nvc, scalestor, offsetstor, &csamp, locstor, offstor, phasorstor, &irow, &dow, &cinv);
  }  

 }//end pragma parallel


 
 timemass_p +=tim.real();
 Int  ixsub, iysub, icounter;
 ixsub=1;
 iysub=1;
  //////@@@@@@@@@@@@@DEBUGGING
  //nth=1;
  ////////@@@@@@@@@@@@@
  if (nth >3){
    ixsub=8;
    iysub=8; 
  }
  else if(nth >1){
     ixsub=2;
     iysub=2; 
  }
  Int rbeg=startRow+1;
  Int rend=endRow+1;
  Block<Matrix<Double> > sumwgt(ixsub*iysub);
  Vector<Double *> swgtptr(ixsub*iysub);
  Vector<Bool> swgtdel(ixsub*iysub);
  for (icounter=0; icounter < ixsub*iysub; ++icounter){
    sumwgt[icounter].resize(sumWeight.shape());
    sumwgt[icounter].set(0.0);
    swgtptr[icounter]=sumwgt[icounter].getStorage(swgtdel(icounter));
  }
  //cerr << "done thread " << doneThreadPartition_p << "  " << ixsub*iysub << endl;
   if(doneThreadPartition_p < 0){
    xsect_p.resize(ixsub*iysub);
    ysect_p.resize(ixsub*iysub);
    nxsect_p.resize(ixsub*iysub);
    nysect_p.resize(ixsub*iysub);
    for (icounter=0; icounter < ixsub*iysub; ++icounter){
      findGridSector(nx, ny, ixsub, iysub, 0, 0, icounter, xsect_p(icounter), ysect_p(icounter), nxsect_p(icounter), nysect_p(icounter), true);
    }
  }
   Vector<Int> xsect, ysect, nxsect, nysect;
   xsect=xsect_p; ysect=ysect_p; nxsect=nxsect_p; nysect=nysect_p;
   //cerr << xsect.shape() << "  " << xsect << endl;
  const Int* pmapstor=polMap.getStorage(del);
  const Int* cmapstor=chanMap.getStorage(del);
// Dummy sumwt for gridweight part
  Matrix<Double> dumSumWeight(npol, nchan);
  dumSumWeight=sumWeight;
  Bool isDSWC;
  Double *dsumwtstor=dumSumWeight.getStorage(isDSWC);
  Int nc=nchan;
  Int np=npol;
  Int nxp=nx;
  Int nyp=ny;
  Int csize=convSize;
  const Int * flagstor=flags.getStorage(del);
  const Int * rowflagstor=rowFlags.getStorage(del);
  const Int *convsupportstor=convSupportPlanes_p.getStorage(del);
  const Int *convrowmapstor=convRowMap_p.getStorage(del);
  const Int *convchanmapstor=convChanMap_p.getStorage(del);
  const Int *convpolmapstor=convPolMap_p.getStorage(del);
  ///

  
  ////////
  tim.mark(); 

  //  if(useDoubleGrid_p) { //always using double prec here 
  {
    DComplex *gridstor=griddedData2.getStorage(gridcopy);
    
#pragma omp parallel default(none) private(icounter, del) firstprivate(idopsf,  datStorage, wgtStorage, flagstor, rowflagstor, convstor, wconvstor, pmapstor, cmapstor, gridstor,  convsupportstor, nxp, nyp, np, nc,ixsub, iysub, rend, rbeg, csamp, csize, nvp, nvc, nvisrow, phasorstor, locstor, offstor, convrowmapstor, convchanmapstor, convpolmapstor, nPolConv, nChanConv, nConvFunc,xsect, ysect, nxsect, nysect) shared(swgtptr) 
    {   
#pragma omp for schedule(dynamic)      
    for(icounter=0; icounter < ixsub*iysub; ++icounter){
      Int x0=xsect(icounter);
      Int y0=ysect(icounter);
      Int nxsub=nxsect(icounter);
      Int nysub=nysect(icounter);
      

    sectgmosd3(datStorage,
	   &nvp,
	   &nvc,
	   &idopsf,
	   flagstor,
	   rowflagstor,
	   wgtStorage,
	   &nvisrow,
	   gridstor,
	   &nxp,
	   &nyp,
	   &np,
	   &nc,
	   convsupportstor, 
	   &csize,
	   &csamp,
	   convstor,
	   cmapstor,
	   pmapstor,
	   swgtptr[icounter],
	   convrowmapstor,
	   convchanmapstor,
	   convpolmapstor,
	       &nConvFunc, &nChanConv, &nPolConv,
	       &x0, &y0, &nxsub, &nysub, &rbeg, &rend, locstor, offstor,
		 phasorstor
	       );
    }
    }//end pragma parallel
    for (icounter=0; icounter < ixsub*iysub; ++icounter){
      sumwgt[icounter].putStorage(swgtptr[icounter],swgtdel[icounter]);
      sumWeight=sumWeight+sumwgt[icounter];
    }    

    //cerr << "SUMWEIG " << sumWeight << endl;
    griddedData2.putStorage(gridstor, gridcopy);
    if(dopsf && (nth >4))
      tweakGridSector(nx, ny, ixsub, iysub);
    timegrid_p+=tim.real();
    tim.mark();
    if(!doneWeightImage_p){
      //This can be parallelized by making copy of the central part of the griddedWeight
      //and adding it after dooing the gridding
      DComplex *gridwgtstor=griddedWeight2.getStorage(weightcopy);
      gmoswgtd2(&nvp, &nvc,flagstor, rowflagstor, wgtStorage, &nvisrow, 
	       &nxp, &nyp, &np, &nc, convsupportstor, &csize, &csamp, 
	       cmapstor, pmapstor,
	       gridwgtstor, dsumwtstor, wconvstor, convrowmapstor, 
	       convchanmapstor,  convpolmapstor, 
	       &nConvFunc, &nChanConv, &nPolConv, &rbeg, 
	       &rend, locstor, offstor, phasorstor);
      griddedWeight2.putStorage(gridwgtstor, weightcopy);
    
    }
    timemass_p+=tim.real();
  }

  convFunc.freeStorage(convstor, convcopy);
  weightConvFunc_p.freeStorage(wconvstor, wconvcopy);
  dumSumWeight.putStorage(dsumwtstor, isDSWC);
  //cerr << "dumSumwe " << dumSumWeight << endl;
  uvw.freeStorage(uvwstor, uvwcopy);
  if(!dopsf)
    data.freeStorage(datStorage, isCopy);

  elWeight.freeStorage(wgtStorage,iswgtCopy);
  



}
*/


/*
void AWPLPG::gridImgWeights(const vi::VisBuffer2& vb){

  if(doneWeightImage_p)
    return;
  matchChannel(vb);
 
  
  //cerr << "CHANMAP " << chanMap << endl;
  //No point in reading data if its not matching in frequency
  if(max(chanMap)==-1)
    return;

  Int startRow, endRow, nRow;
  nRow=vb.nRows();
  startRow=0;
  endRow=nRow-1;
  
  
  //const Matrix<Float> *imagingweight;
  //imagingweight=&(vb.imagingWeight());
  Matrix<Float> imagingweight;
  getImagingWeight(imagingweight, vb);


  Cube<Complex> data;
  //Fortran gridder need the flag as ints 
  Cube<Int> flags;
  Matrix<Float> elWeight;
  interpolateFrequencyTogrid(vb, imagingweight,data, flags, elWeight, FTMachine::PSF);
  
 

  Bool iswgtCopy;
  const Float *wgtStorage;
  wgtStorage=elWeight.getStorage(iswgtCopy);
  Bool issumWgtCopy;
  Double* sumwgtstor=sumWeight.getStorage(issumWgtCopy);

  
 
  // Get the uvws in a form that Fortran can use and do that
  // necessary phase rotation. 
  Matrix<Double> uvw(negateUV(vb));
  Vector<Double> dphase(vb.nRows());
  dphase=0.0;
 
  doUVWRotation_p=true;
  girarUVW(uvw, dphase, vb);
  refocus(uvw, vb.antenna1(), vb.antenna2(), dphase, vb);
  // This needs to be after the interp to get the interpolated channels
  //Also has to be after rotateuvw in case tracking is on
  findConvFunction(*image, vb);
  //nothing to grid here as the pointing resulted in a zero support convfunc
  if(convSupport <= 0)
    return;
  
  Bool del;
  
  const Int* pmapstor=polMap.getStorage(del);
  const Int* cmapstor=chanMap.getStorage(del);
  
  Vector<Int> rowFlags(vb.nRows());
  rowFlags=0;
  rowFlags(vb.flagRow())=true;
  if(!usezero_p) {
    for (uInt rownr=0; rownr< vb.nRows(); rownr++) {
      if(vb.antenna1()(rownr)==vb.antenna2()(rownr)) rowFlags(rownr)=1;
    }
  }

  //Fortran indexing
  
  Int rbeg=1;
  Int rend=vb.nRows();

  const Int * flagstor=flags.getStorage(del);
  const Int * rowflagstor=rowFlags.getStorage(del);

  const Int *convrowmapstor=convRowMap_p.getStorage(del);
  const Int *convchanmapstor=convChanMap_p.getStorage(del);
  const Int *convpolmapstor=convPolMap_p.getStorage(del);
  const Int *convsupportstor=convSupportPlanes_p.getStorage(del);
  //Tell the gridder to grid the weights too ...need to do that once only
  //Int doWeightGridding=1;
  //if(doneWeightImage_p)
  //  doWeightGridding=-1;
  //    IPosition s(flags.shape());
  const IPosition& fs=flags.shape();
  //cerr << "flags shape " << fs << endl;
  std::vector<Int>s(fs.begin(), fs.end());
  Int nvp=s[0];
  Int nvc=s[1];
  Int nvisrow=s[2];
  Int csamp=convSampling;
  Bool uvwcopy; 
  const Double *uvwstor=uvw.getStorage(uvwcopy);
  Bool gridcopy;
  Bool convcopy;
  Bool wconvcopy;
  const Complex *wconvstor=weightConvFunc_p.getStorage(wconvcopy);
  Int nPolConv=convFunc.shape()[2];
  Int nChanConv=convFunc.shape()[3];
  Int nConvFunc=convFunc.shape()(4);
  Bool weightcopy;
  ////////@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  Cube<Int> loc(2, nvc, vb.nRows());
  Cube<Int> off(2, nvc, vb.nRows());
  Matrix<Complex> phasor(nvc, vb.nRows());
  Bool delphase;
  Complex * phasorstor=phasor.getStorage(delphase);
  const Double * visfreqstor=interpVisFreq_p.getStorage(del);
  const Double * scalestor=uvScale.getStorage(del);
  const Double * offsetstor=uvOffset.getStorage(del);
  Int * locstor=loc.getStorage(del);
  Int * offstor=off.getStorage(del);
  const Double *dpstor=dphase.getStorage(del);

  Int irow;
  Int nth=1;
#ifdef _OPENMP
  if(numthreads_p >0){
    nth=min(numthreads_p, omp_get_max_threads());
  }
  else{   
    nth= omp_get_max_threads();
  }
  //nth=min(4,nth);
#endif

  Double cinv=Double(1.0)/C::c;
 
  Int dow=0;

#pragma omp parallel default(none) private(irow) firstprivate(visfreqstor, nvc, scalestor, offsetstor, csamp, phasorstor, uvwstor, locstor, offstor, dpstor, dow, cinv) shared(startRow, endRow) num_threads(nth)  
{
#pragma omp for
  for (irow=startRow; irow<=endRow;irow++){
    
    locuvw(uvwstor, dpstor, visfreqstor, &nvc, scalestor, offsetstor, &csamp, locstor, offstor, phasorstor, &irow, &dow, &cinv);
  }  

 }//end pragma parallel



//always using double prec in this gridder
//  if(useDoubleGrid_p) {
 {
      //This can be parallelized by making copy of the central part of the griddedWeight
      //and adding it after dooing the gridding
      DComplex *gridwgtstor=griddedWeight2.getStorage(weightcopy);
      gmoswgtd2(&nvp, &nvc,flagstor, rowflagstor, wgtStorage, &nvisrow, 
	       &nx, &ny, &npol, &nchan, convsupportstor, &convSize, &convSampling, 
	       cmapstor, pmapstor,
	       gridwgtstor, sumwgtstor, wconvstor, convrowmapstor, 
	       convchanmapstor,  convpolmapstor, 
	       &nConvFunc, &nChanConv, &nPolConv, &rbeg, 
	       &rend, locstor, offstor, phasorstor);
      griddedWeight2.putStorage(gridwgtstor, weightcopy);
    
    



  }

  sumWeight.putStorage(sumwgtstor, issumWgtCopy); 
  elWeight.freeStorage(wgtStorage,iswgtCopy);
    
}
*/
/*
void AWPLPG::get(vi::VisBuffer2& vb, Int row)
{
  

  
  // If row is -1 then we pass through all rows
  Int startRow, endRow, nRow;
  if (row==-1) {
    nRow=vb.nRows();
    startRow=0;
    endRow=nRow-1;
    //  vb.modelVisCube()=Complex(0.0,0.0);
  } else {
    nRow=1;
    startRow=row;
    endRow=row;
    //  vb.modelVisCube().xyPlane(row)=Complex(0.0,0.0);
  }
  

 

  matchChannel(vb);
 
  //No point in reading data if its not matching in frequency
  if(max(chanMap)==-1)
    return;

  // Get the uvws in a form that Fortran can use
  Matrix<Double> uvw(negateUV(vb));
  Vector<Double> dphase(vb.nRows());
  dphase=0.0;
 
  doUVWRotation_p=true;
  girarUVW(uvw, dphase, vb);
  refocus(uvw, vb.antenna1(), vb.antenna2(), dphase, vb);
  
  
  
 
  Cube<Complex> data;
  Cube<Int> flags;
  getInterpolateArrays(vb, data, flags);
  
  //Need to get interpolated freqs
  findConvFunction(*image, vb);

  // no valid pointing in this buffer
  if(convSupport <= 0)
    return;
  Complex *datStorage;
  Bool isCopy;
  datStorage=data.getStorage(isCopy);
  

  Vector<Int> rowFlags(vb.nRows());
  rowFlags=0;
  rowFlags(vb.flagRow())=true;
  if(!usezero_p) {
    for (Int rownr=startRow; rownr<=endRow; rownr++) {
      if(vb.antenna1()(rownr)==vb.antenna2()(rownr)) rowFlags(rownr)=1;
    }
  }
  Int nvp=data.shape()[0];
  Int nvc=data.shape()[1];
  Int nvisrow=data.shape()[2];
  Int csamp=convSampling;
  Int csize=convSize;
  //Int csupp=convSupport;
  Int nc=nchan;
  Int np=npol;
  Int nxp=nx;
  Int nyp=ny;
  Bool uvwcopy; 
  const Double *uvwstor=uvw.getStorage(uvwcopy);
  Int nPolConv=convFunc.shape()[2];
  Int nChanConv=convFunc.shape()[3];
  Int nConvFunc=convFunc.shape()(4);
  ////////@@@@@@@@@
  Cube<Int> loc(2, nvc, nRow);
  Cube<Int> off(2, nvc, nRow);
  Matrix<Complex> phasor(nvc, nRow);
  Bool delphase;
  Bool del;
  const Int* pmapstor=polMap.getStorage(del);
  const Int* cmapstor=chanMap.getStorage(del);
  Complex * phasorstor=phasor.getStorage(delphase);
  const Double * visfreqstor=interpVisFreq_p.getStorage(del);
  const Double * scalestor=uvScale.getStorage(del);
  const Double * offsetstor=uvOffset.getStorage(del);
  const Int * flagstor=flags.getStorage(del);
  const Int * rowflagstor=rowFlags.getStorage(del);
  Int * locstor=loc.getStorage(del);
  Int * offstor=off.getStorage(del);
  const Double *dpstor=dphase.getStorage(del);
  const Int *convrowmapstor=convRowMap_p.getStorage(del);
  const Int *convchanmapstor=convChanMap_p.getStorage(del);
  const Int *convpolmapstor=convPolMap_p.getStorage(del);
  const Int *convsupportstor=convSupportPlanes_p.getStorage(del);
  ////////@@@@@@@@@@@@@@@@@@@@@

  Int irow;
  Int nth=1;
 #ifdef _OPENMP
  if(numthreads_p >0){
    nth=min(numthreads_p, omp_get_max_threads());
  }
  else{   
    nth= omp_get_max_threads();
  }
  //nth=min(4,nth);
#endif
 
  Timer tim;
  tim.mark();

   Int dow=0;
   Double cinv=Double(1.0)/C::c;
#pragma omp parallel default(none) private(irow) firstprivate(visfreqstor, nvc, scalestor, offsetstor, csamp, phasorstor, uvwstor, locstor, offstor, dpstor, dow, cinv) shared(startRow, endRow) num_threads(nth)  
{
#pragma omp for
  for (irow=startRow; irow<=endRow;irow++){
    /////////////////locateuvw(uvwstor,dpstor, visfreqstor, nvc, scalestor, offsetstor, csamp, 
    //    locstor, 
		///////////	      offstor, phasorstor, irow, false);
    //using the fortran version which is significantly faster ...this can account for 10% less overall degridding time
    locuvw(uvwstor, dpstor, visfreqstor, &nvc, scalestor, offsetstor, &csamp, locstor, 
	   offstor, phasorstor, &irow, &dow, &cinv);
  }  

 }//end pragma parallel
 Int rbeg=startRow+1;
 Int rend=endRow+1;
 Int npart=nth;
 
 Bool gridcopy;
 const Complex *gridstor=griddedData.getStorage(gridcopy);
 Bool convcopy;
 ////Degridding needs the conjugate ...doing it here
 Array<Complex> conjConvFunc=conj(convFunc);
 const Complex *convstor=conjConvFunc.getStorage(convcopy);
  Int ix=0;
#pragma omp parallel default(none) private(ix, rbeg, rend) firstprivate(uvwstor, datStorage, flagstor, rowflagstor, convstor, pmapstor, cmapstor, gridstor, nxp, nyp, np, nc, csamp, csize, convsupportstor, nvp, nvc, nvisrow, phasorstor, locstor, offstor, nPolConv, nChanConv, nConvFunc, convrowmapstor, convpolmapstor, convchanmapstor, npart)  num_threads(npart)
  {
    #pragma omp for schedule(dynamic) 
    for (ix=0; ix< npart; ++ix){
      rbeg=ix*(nvisrow/npart)+1;
      rend=(ix != (npart-1)) ? (rbeg+(nvisrow/npart)-1) : (rbeg+(nvisrow/npart)-1+nvisrow%npart) ;
      //cerr << "maps "  << convChanMap_p << "   " << chanMap  << endl;
      //cerr << "nchan " << nchan << "  nchanconv " << nChanConv << " npolconv " << nPolConv << " nRowConv " << nConvFunc << endl;
     sectdmos3(
	       datStorage,
	       &nvp,
	       &nvc,
	       flagstor,
	       rowflagstor,
	       &nvisrow,
	       gridstor,
	       &nxp,
	       &nyp,
	       &np,
	       &nc,
	       convsupportstor,
	       &csize,   
	       &csamp,
	       convstor,
	       cmapstor,
	       pmapstor,
	       convrowmapstor, convchanmapstor,
	       convpolmapstor,
	       &nConvFunc, &nChanConv, &nPolConv,
	       &rbeg, &rend, locstor, offstor, phasorstor
	       );


    }
  }//end pragma omp


  data.putStorage(datStorage, isCopy);
  griddedData.freeStorage(gridstor, gridcopy);
  convFunc.freeStorage(convstor, convcopy);
  
   timedegrid_p+=tim.real();

  interpolateFrequencyFromgrid(vb, data, FTMachine::MODEL);
}
*/

  } // REFIM ends
} //# NAMESPACE CASA - END
