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

  
  AWPLPG::AWPLPG(SkyJones* sj, const Int nw,  Bool dosquint, const Double painc, MPosition mloc, String stokes,  const Bool usezero, const Bool useDoublePrec,  const casacore::Bool usePointing): MosaicFTNew(sj, mloc, stokes, Long(1000000), 16, usezero, useDoublePrec, False, usePointing), doSquint_p(dosquint), paInc_p(painc), nw_p(nw) {



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
              pbConvFunc_p->findUsefulChannels(chunkfreq, vb);
              cerr <<  "chunkfreq " <<  chunkfreq <<  endl;
              std::move(chunkfreq.begin(), chunkfreq.end(), std::back_inserter(freqs));
              if(doSquint_p)
                pAs.push_back(getPA(vb));
              if(nw_p > 1)
                	maxW=max(maxW, max(abs(vb.uvw().row(2)*max(vb.getFrequencies(0))))/C::c);
          }
    }
    
    //return vi to origin
    vi->originChunks(); vi->origin();
    std::sort(freqs.begin(),  freqs.end());
    auto last = std::unique(freqs.begin(),  freqs.end());
    freqs.erase(last,  freqs.end());
    Double paInc=0.0;
    if(pAs.size()>1){
      std::sort(pAs.begin(), pAs.end());
      last=std::unique(pAs.begin(), pAs.end());
      pAs.erase(last, pAs.end());
      if(pAs.size() >1)
        paInc=fabs(pAs[0]-pAs[1]);
    }
    
    cerr <<  "Freqs " <<  freqs <<  endl;
    if (nw_p == 0)
      nw_p = 1;
    Vector<Double> wVals(nw_p,0);
    if(nw_p >1){
      Double st=maxW/(Double(nw_p-1)*Double(nw_p-1));
      for (int k=0; k <nw_p; ++k)
        wVals[k]=Double(k*k)*st;
    }
    (*awConvs_p).addConvFunc(Vector<Double>(freqs), wVals, paInc);
    
  }
  
}
  
 void AWPLPG::findConvFunction(const ImageInterface<Complex>& iimage, const vi::VisBuffer2& vb ){
  //
  // pbConvFunc_p.phasegradient
    convFunc.resize();
    convFunc.assign(awConvs_p->getConvFunc());
    /* { 
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
      PagedImage<Complex> lastplane(elshp,  csysA,  "NOOBOO");
      lastplane.put(convFunc(elblc,  eltrc).nonDegenerate());
    
    //////
    } */  
    weightConvFunc_p.resize();
    weightConvFunc_p.assign(awConvs_p->getWeightConvFunc());
    convSizePlanes_p.resize();
    convSizePlanes_p = awConvs_p->getConvSizes();
    convSupportPlanes_p.resize();
    convSupportPlanes_p = awConvs_p->getConvSupports();
    awConvs_p->getConvIndices(convPolMap_p,  convChanMap_p,  convRowMap_p,  vb);
    //cerr <<  "min max convrowmap " <<  min(convRowMap_p) <<  "  " <<  max(convRowMap_p) <<  " supp " <<   max(convSupportPlanes_p) <<  " csize " << max(convSizePlanes_p) <<  " convchanmap "<< min(convChanMap_p) <<  "    " << max(convChanMap_p) << endl;
    
    pbConvFunc_p->rephaseConvFunc(iimage, vb, convSampling,  convFunc, weightConvFunc_p,  MVDirection(-(movingDirShift_p.getAngle())), fixMovingSource_p);
    convSupport =max(convSupportPlanes_p);
    convSize = max(convSizePlanes_p);
   
 }


  } // REFIM ends
} //# NAMESPACE CASA - END
