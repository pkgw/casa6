
// -*- C++ -*-
//# AWVisResampler.cc: Implementation of the AWVisResampler class
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
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$
#ifdef USE_HPG
#include <synthesis/TransformMachines/SynthesisError.h>
#include <synthesis/TransformMachines2/AWVisResamplerHPG.h>
#include <synthesis/TransformMachines2/AWConvFuncHolder.h>
#include <synthesis/TransformMachines2/Utils.h>
#include <synthesis/TransformMachines/SynthesisMath.h>
#include <casacore/coordinates/Coordinates/SpectralCoordinate.h>
#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/casa/OS/Timer.h>
#include <fstream>
#include <iostream>
#include <typeinfo>
#include <iomanip>
#include <map>
//#include <synthesis/TransformMachines2/FortranizedLoops.h>
//#include <synthesis/TransformMachines2/hpg.hpp>

//#include <hpg/hpg.hpp>
//#include <hpg/hpg_indexing.hpp>

#include <type_traits>
#include <synthesis/TransformMachines2/HPGVisBuffer.cc>

#ifdef _OPENMP
#include <omp.h>
#endif

//#include <casa/BasicMath/Functors.h>
using namespace casacore;
using namespace hpg;
namespace casa{
  using namespace refim;

  //
  // Global functions
  //
  // ---------------
  
  template <unsigned N> std::vector<hpg::VisData<N>>
makeHPGVisBuffer2(casacore::Matrix<double>& sumwt,
		 const casa::refim::VBStore& casaVBS,
		 AWConvFuncHolder& awh,
		 const unsigned nGridPol, const unsigned nGridChan,
		 unsigned startRow,  unsigned endRow,
		 const unsigned startChan, const unsigned endChan,
		 const casacore::Vector<casacore::Int>& chanMap,
		 const casacore::Vector<casacore::Int>& polMap,
		 const casacore::Vector<double>& dphase,
		 const hpg::CFSimpleIndexer& cfsi
		 )
{
	unsigned nVisChan=endChan - startChan + 1;
	unsigned nVisRow=endRow - startRow + 1;
	const casa::VisBuffer2& vb = *(casaVBS.vb_p);
	IPosition dataShape = vb.visCube().shape();
	unsigned nDataPol=dataShape(0);
	vector<int> unique_pol = polMap.tovector();
	std::sort(unique_pol.begin(),  unique_pol.end());
	auto last = std::unique(unique_pol.begin(),  unique_pol.end());
	unique_pol.erase(last,  unique_pol.end());
	uint nVisPol = unique_pol.size();
	
	
	std::vector<hpg::VisData<N>> hpgVB=create_blank_vis_data_vector<N>(nVisPol,nVisChan,nVisRow);
	std::array<std::complex<hpg::visibility_fp>, N> vis;
	std::array<hpg::vis_weight_fp, N> wt;
	
	const casacore::Matrix<double> UVW=casaVBS.uvw_p;
	
	std::array<hpg::vis_uvw_fp, 3> hpgUVW;
	hpg::vis_phase_fp              dphaseHPG;
	unsigned                       grid_cube=0; // For now, grid on the same plane.
	// Need to pass the beam offsets phase grad
	//setting to use pointing table for now
	Vector<Double> pointingOffsets = awh.getPointingPhaseShift(vb, True);
	hpg::cf_phase_gradient_t       cf_phase_gradient = {(hpg::cf_phase_gradient_fp)pointingOffsets[0],
					 (hpg::cf_phase_gradient_fp)pointingOffsets[1]};
	   
	
	std::array<unsigned, 2>        cf_index;
	Vector<Double>freq = vb.getFrequencies(0);
//	Int vbSpw = vb.spectralWindows()(0);
	unsigned hpgIRow=0;
	// Lets get convindices
	Vector<Int> convpolmap;
	Vector<Int> convchanmap;
	Vector<Int> convrowmap;
	awh.getConvIndicesHPG(convpolmap,  convchanmap,  convrowmap,  vb,  UVW);
	// convrowmap goes from 0 to nW-1 inclusive
	uint nconvchan = awh.getFreqValsHPG().nelements();
	uint nconvwvals = awh.getWValsHPG().nelements();
	//cerr <<  "vb.spec" << vb.spectralWindows()(0) <<  "convchanmap" <<  convchanmap <<  " nconvchan " << nconvchan << " convrow " << convrowmap <<   endl;
	for(unsigned irow=startRow; irow< endRow; irow++) 
	 {
	  if (casaVBS.ftmType_p==casa::refim::FTMachine::WEIGHT)
		{
			  hpgUVW={0.0, 0.0, 0.0};
		}
		else{
		  hpgUVW={UVW(0,irow),
				  UVW(1,irow),
				  UVW(2,irow)};
		}
		
	   if(!casaVBS.rowFlag_p(irow) )
		{
		 for(unsigned ichan=startChan; ichan< endChan; ichan++)
	    {
		 if(((chanMap[ichan]>=0) && (chanMap[ichan]<nGridChan))) 
			{
				int polelem=0;  
			for(unsigned ipol=0; ipol< nDataPol; ipol++) 
			{
			    
				if ((polMap(ipol)>=0) && (polMap(ipol)<nGridPol)) 
				{
				      //int visVecElement= polMap[ipol];
					  int visVecElement=polelem;
				      if (vb.flagCube()(ipol,ichan,irow)==false)
					{
					  if ((casaVBS.ftmType_p==casa::refim::FTMachine::WEIGHT)||
					      (casaVBS.ftmType_p==casa::refim::FTMachine::PSF))
					    {
					      vis[visVecElement] = std::complex<double>(1.0,0.0);
					    }
					  else
					    {
					      vis[visVecElement] = casaVBS.visCube_p(ipol,ichan,irow);
					    }

					  wt[visVecElement]=casaVBS.imagingWeight_p(ichan,irow);
					}
				      else
					  {
					  vis[visVecElement] = 0.0;
					  wt[visVecElement]  = 0.0;
					}
					++polelem;
				  	sumwt(polMap[ipol],chanMap[ichan]) += vb.imagingWeight()(ichan, irow);
				}
				// if polmap
			    
			}
			//ipol
			uint cindex=0;
			hpg::vis_frequency_fp frequency = freq[ichan];
			if (casaVBS.ftmType_p==casa::refim::FTMachine::PSF || casaVBS.ftmType_p==casa::refim::FTMachine::WEIGHT ) 
			{
				dphaseHPG = 0.0;
				//cindex = (0 + nconvwvals - 1) * nconvchan + convchanmap[ichan];

			}
			else{
				dphaseHPG = -2.0 * C::pi * dphase[irow] * frequency / C::c;
				//cindex = (convrowmap[irow] + nconvwvals - 1) * nconvchan + convchanmap[ichan];
			}
			//uint cindex = (convrowmap[irow]+nconvwvals-1)*nconvchan+convchanmap[ichan];
			cindex = convchanmap[ichan] * (nconvwvals) + abs(convrowmap[irow]);

			cf_index = {0,cindex};
			// cfindex is rowindex*nconvchan+ convchanmap
			hpgVB[hpgIRow++] = hpg::VisData<N>(vis,wt,frequency,dphaseHPG,hpgUVW,grid_cube,cf_index
							 ,cf_phase_gradient
							 );
			}
		// if chnamap
		}
		// ichan
  }													                                            // if flag
	}                                                          // irow
  hpgVB.resize(hpgIRow);
  return hpgVB;
}

  
  
  
  
  
  
  
 //----------------------------------------------------------------------------------//----
  template <unsigned N>
  void makeMuellerIndexes(const PolMapType& mNdx,
			  const PolMapType& mValues,
			  const unsigned& nGridPol,
			  std::vector<std::array<int, N>>& muellerIndexes)
  {
    PolMapType tt;
    //cerr << "####N " << N << endl;
    // for(unsigned i=0;i<mNdx.size();i++)
    //   {
    // 	cerr << mNdx[i] << endl;
    //   }
    // cerr << endl;

    for(unsigned i=0,k=0;i<mNdx.size();i++)
      {
	if (mNdx[i].size()>0)
	  {
	    tt.resize(k+1,true);
	    tt[k].resize(mNdx[i].size());
	    tt[k]=mNdx[i];
	    k++;
	  }
      }

    muellerIndexes.resize(nGridPol);
    for(unsigned i=0;i<muellerIndexes.size();i++)
      for(unsigned j=0;j<muellerIndexes[i].size();j++)
	muellerIndexes[i][j]=-1;
    ////This is incomprehensible indexing (2023/3/13) but mValues size and tt size may not match ...so for now looping
    //// over the minimum of the two sizes just to avoid out of bounds indexing !!
    for(unsigned i=0 ; i<min(mValues.size(), tt.size()) ;i++)
      {
        //cerr <<"i= " << i <<  " SIZES mValues " << mValues[i].size() << " mval[i] "<< mValues[i] << " tt[i] " << tt[i] << endl;
	for(unsigned j=0;j<mValues[i].size();j++)
	  {
	    muellerIndexes[j][mValues[i][j]%N]=tt[i][j];
	    //muellerIndexes[j][mValues[i][j]%N]=1;
	  }
      }
    /*
    for(unsigned ir=0;ir<muellerIndexes.size();ir++)
       {
     	for(unsigned ic=0;ic<muellerIndexes[ir].size();ic++)
     	    cerr << muellerIndexes[ir][ic] << " ";
     	cerr << endl;
       }
    */
  }
  //
  //-------------------------------------------------------------------------
  //  
  void makeAWLists(CFBuffer& cfb, const bool& wbAWP, const int& nw,
		   const double& imRefFreq, const double& spwRefFreq,
		   Vector<Int>& wNdxList, Vector<Int>& spwNdxList,
		   const int vbSPW)
  {
    Vector<Double> wVals, fVals; PolMapType mVals, mNdx, conjMVals, conjMNdx;
    Double fIncr, wIncr;

    //
    // The following can be generalized to pick a subset of CFs along
    // W and SPW axis in the CFB.  Perhaps also useful in the longer
    // run, e.g. when using a super-set CFC not all of which may be
    // used for a given imaging.
    //
    cfb.getCoordList(fVals,wVals,mNdx, mVals, conjMNdx, conjMVals, fIncr, wIncr);

    // W-pixels in the CFC should be >= w-planes user setting
    assert(wVals.nelements() >= nw);

    // Make list of W-CF indexes
    int nWCFs=(nw==1)?1:wVals.nelements();
    wNdxList.resize(nWCFs);
    for(int i=0;i<nWCFs;i++) wNdxList[i] = i;


    //   cerr << "nWCFs" << nWCFs << " vbSPW " << vbSPW << " fvals " << fVals << endl;
    // Make list of SPW-CF indexes
    int nSPWCFs=fVals.nelements();
    if (wbAWP==true)
      {
	// If a valid SPW ID is given, translate it to the spwNdx for the nearest SPW
	if ((vbSPW>=0))// && (vbSPW <nSPWCFs))
	  {
	    int refSPW;
	    std::vector<double> stdList(nSPWCFs);
	    for (int i=0; i<nSPWCFs; i++) stdList[i] = fVals[i];
	    Double refCFFreq = SynthesisUtils::stdNearestValue(stdList, spwRefFreq, refSPW);

	    spwNdxList.resize(1);
	    spwNdxList[0]=refSPW;
	  }
	else
	  {
	    spwNdxList.resize(nSPWCFs);
	    for(int i=0;i<nSPWCFs;i++) spwNdxList[i] = i;
	  }
      }
    else
      {
	// For wbAWP=F, pick up the CF closest to the image reference frequency
	int refSPW;
	std::vector<double> stdList(nSPWCFs);
	for (int i=0; i<nSPWCFs; i++) stdList[i] = fVals[i];
	Double refCFFreq = SynthesisUtils::stdNearestValue(stdList, imRefFreq, refSPW);

	spwNdxList.resize(1);
	spwNdxList[0]=refSPW;
       
      }
    //  cerr << "spwNdxList " << spwNdxList << endl;
    return;
  }
  //
  //-------------------------------------------------------------------------
  std::tuple<hpg::opt_t<hpg::Error>, hpg::CFSimpleIndexer> 
    AWVisResamplerHPG::loadCF(VBStore& vbs, AWConvFuncHolder& awh, 
			  bool send_to_device) {
      
      LogIO log_l(LogOrigin("AWVisResamplerHPG","loadCF"));
      Int ndatapol = vbs.vb().nCorrelations();
      Int ndataChan = vbs.vb().nChannels();
      awh.resetHPGConvFuncs(vbs.vb());
      Array<Complex> cFunc;
      if (vbs.ftmType_p !=  casa::refim::FTMachine::WEIGHT) {
		cFunc = awh.getConvFuncHPG();
	  }
	   else{		
		cFunc = awh.getWeightConvFuncHPG();
	   }
      IPosition cshape = cFunc.shape();
      cerr <<  "convfunc shape" <<  cshape <<  " npix " <<  cshape.product() <<  endl;
      uint convSize = cshape[0];
      uInt nConvPol = cshape[2];
      Vector<Double> cFreqs = awh.getFreqValsHPG();
      Vector<Double> cWVals = awh.getWValsHPG();
      uint nPol = 2;                                         // we'll try only the RR
      uint nchan = cFreqs.nelements();
      uint nw = cWVals.nelements();
      int sampling = awh.getOverSampling();
      uint nCF = nchan * (nw);
       //                               nBLType, nTime/PA, nW, nFreq, nPol
      hpg::CFSimpleIndexer cfsi({1,false},{1,false},{nw,true},{nchan,true}, nPol);

      if (cfArray.oversampling()==0) {
        
        log_l << "Setting cfArray size: " << "nCF: " << nCF << " x " << nPol << " sampling: " << sampling <<  " convsize " <<  convSize << LogIO::POST;
        cfArray.setSize((unsigned)(nCF),(unsigned)sampling);
        
      }
      Bool isCopy;
      Complex *cfuncPtr = cFunc.getStorage(isCopy);
      //Matrix<Complex> conjW(convSize,  convSize);
      //Bool isConjCopy;
	  //Complex *conjptr = conjW.getStorage(isConjCopy);
      unsigned iGrp=0; // HPG Group index
      cerr <<  "@@@@@LoadCF nchan " <<  nchan <<  " nw " <<  nw <<  endl;
      for(uint iFreq=0; iFreq < nchan; ++iFreq) // CASA CF Freq-index
      {
		  // CASA CF W-index
		  uint wcounter = 0;
		  for (Int iW =0; iW < Int(nw); ++iW)
        {
		  
          for(uint ipol=0; ipol < nPol; ++ipol)  // CASA CF Pol-index
          {
            Complex* convptr=NULL;
			uint ptroffset =  ((abs(iW)*nchan+iFreq)*nConvPol+ipol)*convSize*convSize;
			//uint ptroffset = ((iFreq * nw + iW) * nConvPol + ipol) * convSize * convSize;
			//cerr <<  "ptroffset " <<  ptroffset <<  endl;
            convptr = (cfuncPtr+ptroffset);
            /*if (iW < 0) {
				
				memcpy(conjptr, convptr, sizeof(Complex)*convSize*convSize);
				conjW.putStorage(conjptr, isConjCopy);
				//conjW = conj(conjW);
				conjptr = conjW.getStorage(isConjCopy);
				convptr = conjptr;
			}*/
            hpg::CFCellIndex cfCellidx(0,0, 0,iFreq,ipol);
			std::array<unsigned, 3> index = cfsi.cf_index(cfCellidx);
			cfArray.resize(iGrp, convSize, convSize);
            if (send_to_device) {
			//cerr << "send iFreq " << iFreq <<  " igrp " <<  iGrp << " indices " <<  index[0] <<  "," <<  index[1] <<  "," <<  index[2] << endl;  
             cfArray.setValues(convptr,  iGrp,  convSize,  convSize,  ipol,  0);
            }
          }                                                 // pol
		  ++iGrp;
		  ++wcounter;
		}// iW
	   
      }                                                     // ifreq
      cFunc.putStorage(cfuncPtr,  isCopy);
      hpg::opt_t<hpg::Error> err;
      if (send_to_device) {
        	err=hpgGridder_p->set_convolution_function(hpg::Device::OpenMP, std::move(cfArray));
      }
      std::tuple<hpg::opt_t<hpg::Error>, hpg::CFSimpleIndexer> retval = std::make_tuple(err, cfsi);
      return retval;
      
      
        
  }
      
      
      
      
  //--------------------------------------------------------------------------
  //  
  // This is a global method, used in AWVRHPG::DataToGrid_impl().
  template <unsigned N>
  hpg::Gridder* AWVisResamplerHPG::initGridder2(const hpg::Device HPGDevice_l,
			     const uInt& NProcs,
			     const hpg::CFArray* cfArray_ptr,
			     const std::array<unsigned, 4>& grid_size,
			     const std::array<double, 2>& grid_scale,
			     const PolMapType& mNdx,
			     const PolMapType& mVals,
			     const PolMapType& conjMNdx,
			     const PolMapType& conjMVals,
			     const int& nAntenna,
			     const int& nChannel
			     // std::vector<std::array<int, N> > mueller_indexes,
			     // std::vector<std::array<int, N> > conjugate_mueller_indexes
			     )
    {
      //
      // FYI:  gridSize{(uInt)nx,(uInt)ny,(uInt)nGridPol,(uInt)nGridChan};
      //
      LogIO log_l(LogOrigin("AWVRHPG", "initGridder2[R&D]"));

      log_l << "grid_size: " << grid_size[0] << " " << grid_size[1] << " "  << grid_size[2] << " "  << grid_size[3] << " " << LogIO::POST;
      log_l << "Per VB:  nAnt: " << nAntenna << ", nChannel: " << nChannel << LogIO::POST;

      if (!hpg::is_initialized()) hpg::initialize();


	{
	  std::vector<std::array<int, N> > mueller_indexes,conjugate_mueller_indexes;

	  makeMuellerIndexes<N>(mNdx, mVals, grid_size[2]/*nGridPol*/, mueller_indexes);
	  makeMuellerIndexes<N>(conjMNdx, conjMVals, grid_size[2]/*nGridPol*/, conjugate_mueller_indexes);


	  //size_t max_visibilities_batch_size = 351*2*1000;
	  size_t max_visibilities_batch_size = (nAntenna*(nAntenna-1)/2)*2*nChannel;
	  hpg::rval_t<hpg::Gridder> g;

	  cerr << "M: " << endl;
	  for(unsigned ir=0;ir<mueller_indexes.size();ir++)
	    {
	      for(unsigned ic=0;ic<mueller_indexes[ir].size();ic++)
		{	
		  //		  mueller_indexes[ir][ic] = 1;
		  cerr << "ir ic " <<  ir <<  " " <<  ic <<  "  " << mueller_indexes[ir][ic] << " ";
		}
	      cerr << endl;
	    }
	  cerr << "M*: " << endl;
	  for(unsigned ir=0;ir<conjugate_mueller_indexes.size();ir++)
	    {
	      for(unsigned ic=0;ic<conjugate_mueller_indexes[ir].size();ic++)
		{
		  //		  conjugate_mueller_indexes[ir][ic] = 0;
		  cerr <<  "ir ic " <<  ir <<  " " <<  ic <<  "  " << conjugate_mueller_indexes[ir][ic] << " ";
		}
	      cerr << endl;
	    }

    // 	  mueller_indexes.resize(N);  conjugate_mueller_indexes.resize(N);

	  	  cerr << "Mueller index shape: " << mueller_indexes.size() << "X" << HPGNPOL << endl;
	  //	  if (hpgDevice=="serial") Device=hpg::Device::Serial;
	  g = hpg::Gridder::create<N>(HPGDevice_l, NProcs, max_visibilities_batch_size,
				 cfArray_ptr, grid_size, grid_scale, mueller_indexes,
				 conjugate_mueller_indexes);
	  //hpgGridder_p = new hpg::Gridder(hpg::get_value(std::move(g)));

	  return new hpg::Gridder(hpg::get_value(std::move(g)));
	}
    };
    
   // This is a global method, used in AWVRHPG::DataToGrid_impl().
  template <unsigned N>
  hpg::Gridder* AWVisResamplerHPG::initGridder3(const hpg::Device HPGDevice_l,
			     const uInt& NProcs,
			     const hpg::CFArray* cfArray_ptr,
			     const std::array<unsigned, 4>& grid_size,
			     const std::array<double, 2>& grid_scale,
			     const int& nAntenna,
			     const int& nChannel
			     )
    {
      //
      // FYI:  gridSize{(uInt)nx,(uInt)ny,(uInt)nGridPol,(uInt)nGridChan};
      //
      LogIO log_l(LogOrigin("AWVRHPG", "initGridder3"));

      log_l << "grid_size: " << grid_size[0] << " " << grid_size[1] << " "  << grid_size[2] << " "  << grid_size[3] << " " << LogIO::POST;
      log_l << "Per VB:  nAnt: " << nAntenna << ", nChannel: " << nChannel << LogIO::POST;

      if (!hpg::is_initialized()) hpg::initialize();


	{
	  

	  //size_t max_visibilities_batch_size = 351*2*1000;
	  size_t max_visibilities_batch_size = (nAntenna*(nAntenna-1)/2)*2*nChannel;
	  hpg::rval_t<hpg::Gridder> g;

	  
	  std::vector<std::array<int, N> > mueller_indexes,  mueller_conj_indexes;
	  if (N == 2) {
			mueller_indexes.push_back({0, 1});
			mueller_conj_indexes.push_back({1, 0});
	   
		}
	  else{
		throw(AipsError("HPG does not deal with "+String::toString(N) + " pol "));
	}
   
	  g = hpg::Gridder::create<N>(HPGDevice_l, NProcs, max_visibilities_batch_size,
				 cfArray_ptr, grid_size, grid_scale, mueller_indexes,
				 mueller_indexes);
	  //hpgGridder_p = new hpg::Gridder(hpg::get_value(std::move(g)));

	  return new hpg::Gridder(hpg::get_value(std::move(g)));
	}
    };  
  //#include "HPGLoadCF.cc"
#include "HPGLoadCFNew.inc"
  //
  // END -- Global functions
  // --------------------------------------------------------------------------------------

  template
  void AWVisResamplerHPG::DataToGridImpl_p(Array<DComplex>& grid, VBStore& vbs, 
  					Matrix<Double>& sumwt,const Bool& dopsf,
  					Bool useConjFreqCF); // __restrict__;
  template
  void AWVisResamplerHPG::DataToGridImpl_p(Array<Complex>& grid, VBStore& vbs, 
					Matrix<Double>& sumwt,const Bool& dopsf,
					Bool useConjFreqCF); // __restrict__;




  void AWVisResamplerHPG::copy(const AWVisResamplerHPG& other)
    {
      AWVisResampler::copy(other);
      vb2CFBMap_p=new VB2CFBMap(other.getVBRow2CFBMap());
    }


  // Moved the accumulateFromGrid() method to file to play with
  // multi-threading it to not clutter this file.  Both versions
  // (threaded and non-threaded) are in this file.
  //#include "accumulateFromGrid.cc"

  std::tuple<String, hpg::Device> AWVisResamplerHPG::getHPGDevice()
  {
    String hpgDevice="cuda";
    hpg::Device Device=hpg::Device::Cuda;
    //hpgDevice=refim::SynthesisUtils::getenv("HPGDEVICE",hpgDevice);
    if (hpgDevice=="serial") Device=hpg::Device::Serial;

    return std::make_tuple(hpgDevice,Device);
  }
  //
  // Gather the grid and sum-of-weights from the Device (external to the memory address space), and shut the device.
  //
  void AWVisResamplerHPG::finalizeToSky(casacore::Array<casacore::DComplex>& griddedData, casacore::Matrix<casacore::Double>& sumwt)
  {

    if(!hpgGridder_p){
      cerr << "Calling ResamplerHPG finalize with null  pointer...just returning " << endl;
      return;
    }
    LogIO log_l(LogOrigin("AWVisResamplerHPG[R&D]","finalizeToSky(DCompelx)"));
    log_l << "Finalizing gridding..." << LogIO::POST;


    hpgGridder_p->shift_grid(ShiftDirection::FORWARD);
    hpgGridder_p->apply_grid_fft();
    hpgGridder_p->shift_grid(ShiftDirection::BACKWARD);

    GatherGrids(griddedData,sumwt);
	//	size_t size;
	 //auto wgtptr = getSumWeightsPtr(size);
   //cerr << "@@@SIZE " << size << " sumwt " << ((size >0) ? String::toString(wgtptr.get()[size - 1]) : "null") << endl;

    //      cerr << endl << "HPG::resetGridder()" << endl;
    hpgGridder_p->reset_grid();
    if (hpgGridder_p != NULL) delete hpgGridder_p;
    hpgGridder_p=NULL;


    if (isHPGCustodian_p) hpg::finalize();

    log_l << "...done finalizing gridding" << LogIO::POST;
  } 
  ////====================================================================
  void AWVisResamplerHPG::GatherGrids(casacore::Array<casacore::DComplex>& griddedData, casacore::Matrix<casacore::Double>& sumwt) 
  {
      LogIO log_l(LogOrigin("AWVisResamplerHPG[R&D]","GatherGrids(DCompelx)"));
      {
      	auto gg=hpgGridder_p->grid_values();
	std::unique_ptr<hpg::GridWeightArray> wts=hpgGridder_p->grid_weights();
	log_l << "Got the grid from the device " << griddedData.shape() << LogIO::POST;
	cerr << std::setprecision(20) << "Weights shape " << wts->extent(0) << " x " << wts->extent(1) << endl;

	hpgSoW_p.resize(wts->extent(0));
	for (unsigned i=0;i<hpgSoW_p.size();i++)
	  hpgSoW_p[i].resize(wts->extent(1));

	for (unsigned i=0;i<hpgSoW_p.size();i++)
	  {
	  for (unsigned j=0;j<hpgSoW_p[i].size();j++)
	    {
	      hpgSoW_p[i][j]=(*wts)(i,j);
	      cerr << (*wts)(i,j) << " ";
	    }
	  cerr << endl;
	  }
	cerr << "Before: hpgSoW_p, sumwt: " << sumwt << endl;
	// hpgSoW_p is of type sumofweight_fp (vector<vector<double>>).  sumWeight is a Matrix.
	// sow[i][*] is a Mueller row. i is the index for the polarization product.
	// wt below is a sum of the weights of all the Mueller elements in a row.
	if (hpgSoW_p.size() > 0)
	  {
	    for (int i=0; i<sumwt.shape()(0); i++)
	      {
		double wt=0.0;
		for (unsigned j=0; j<hpgSoW_p[i].size(); j++)
		  wt += hpgSoW_p[i][j];
		for (unsigned j=0; j<sumwt.shape()(1); j++)
		  {
		    sumwt(i,j)=wt;
		  }
	      }
	  }
	cerr << "After: hpgSoW_p, sumwt: " << sumwt << endl;
      	// for (size_t i = 0; i < 4; ++i)
      	//   cerr << gg->extent(i) << " ";
      	IPosition shp = griddedData.shape(),ndx(4);
	Bool dummy;
	casacore::DComplex *stor = griddedData.getStorage(dummy);
	DComplex val, max=casacore::DComplex(0.0);

	// std::memcpy((void *)(griddedData.getStorage(dummy)),
	// 	    (void *)(gg->data()),
	// 	    shp.product()*sizeof(casacore::DComplex));

        // To TEST: All the following nested loops could be replaced with
        // gg->copy_to(hpg::Device::OpenMP, stor);
	unsigned long k=0;
	for(ndx[3]=0;  ndx[3]<shp[3];  ndx[3]++)
	  for(ndx[2]=0;  ndx[2]<shp[2];  ndx[2]++)
	    for(ndx[1]=0;  ndx[1]<shp[1];  ndx[1]++)
	      for(ndx[0]=0;  ndx[0]<shp[0];  ndx[0]++)
		{
		  val = (*gg)(ndx[0],ndx[1],ndx[2],ndx[3]);
		  stor[k++] = val;
		}
	log_l << "Copied grid to CPU buffer " << dummy << " " << max << LogIO::POST;
      }
  };
  //
  // This is the method that can be used via the imaging framework.  A
  // dummy for now.  Perhaps the global function below is a better,
  // design-wise. In which case it should move to SynthesisUltils, so
  // that it can be used wherever necessary in the system.
  //
  void AWVisResamplerHPG::initGridder(const uInt&,//NProcs,
				      const hpg::CFArray*,// cfArray_ptr,
				      const std::array<unsigned, 4>&,//grid_size,
				      const std::array<float, 2>&// grid_scale
				      )
  {}


  void AWVisResamplerHPG::setModelImage(std::shared_ptr<casacore::ImageInterface<casacore::Complex> > modelim){
    modelImage_p=modelim;
  }
  bool AWVisResamplerHPG::sendModelImage(
					 const std::array<unsigned, 4>& gridSize)
  {
    LogIO log_l(LogOrigin("AWVisResamplerHPG[R&D]","sendModelImage"));

    if (HPGModelImageName_p == "" && !modelImage_p)
      {
	//log_l << "Name for model image is blank!  No model image sent to the device" << LogIO::WARN << LogIO::POST;
	return false;
      }
      
    // Open the CASA image and transfer it to hpg::GridValueArray in
    // a code-block so that the CASA image storage is released at the end
    // of the block.
    casacore::CoordinateSystem csys;
    casacore::TempImage<casacore::DComplex> modelImageGrid;
    if(Table::isReadable(HPGModelImageName_p)){
      log_l << "Loading model image " << "\"" << HPGModelImageName_p << "\"" <<  LogIO::POST;

      csys=SynthesisUtils::makeModelGridFromImage(HPGModelImageName_p, modelImageGrid);
    }
    else{
      csys=modelImage_p->coordinates();
      modelImageGrid.resize(modelImage_p->shape());
      modelImageGrid.setCoordinateInfo(csys);
      modelImageGrid.copyData((LatticeExpr<DComplex>) (*modelImage_p));
      
    }
      // Stop if the model image sname is not the same as the grid shape
      {
	for (auto i : {0,1,2,3})
	  if (modelImageGrid.shape()[i] != gridSize[i])
	    {
	      log_l << "Model image shape: " << modelImageGrid.shape() << ". "
		    << "Grid shape: [" << gridSize[0] << "," << gridSize[1] << "," << gridSize[2] << "," << gridSize[3] << "]"
		    << LogIO::WARN;
	      log_l << "Model image is of incorrect shape " << LogIO::EXCEPTION << LogIO::POST;
	    }
      }
      Bool dummy;
      assert(hpg::GridValueArray::rank == modelImageGrid.shape().nelements());
      std::array<unsigned,hpg::GridValueArray::rank> extents={(unsigned)modelImageGrid.shape()[0],
							      (unsigned)modelImageGrid.shape()[1],
							      (unsigned)modelImageGrid.shape()[2],
							      (unsigned)modelImageGrid.shape()[3]};
      casacore::Array<casacore::DComplex> tarr=modelImageGrid.get();
      casacore::DComplex *tstor = tarr.getStorage(dummy);

      std::unique_ptr<hpg::GridValueArray> HPGModelImage;

      HPGModelImage = hpg::GridValueArray::copy_from("whatever",
						     HPGDevice_p, // target_device
						     //hpg::Device::Cuda, //target_device
						     hpg::Device::OpenMP,//Device host_device,
						     (hpg::GridValueArray::value_type*)tstor,
						     extents);


    unsigned shape[hpg::GridValueArray::rank];
      
    log_l << "Sending model image of shape ";
    for (unsigned i=0;i<hpg::GridValueArray::rank;i++)
      {
	shape[i]=HPGModelImage->extent(i);
	log_l << shape[i] << " ";
      }
    log_l << LogIO::POST;

    hpg::opt_t<hpg::Error>
      err = hpgGridder_p->set_model(hpg::Device::OpenMP,std::move(*HPGModelImage));
    //    bool success=(err == nullptr);
    if (err)
    	log_l << "Failed hpg::Gridder::set_model(). Error type: " << static_cast<int>(err->type()) << LogIO::SEVERE;
    else
      log_l << "Finished sending image to the device" << LogIO::POST;
	

    grid_value_fp norm = 1.0;//((grid_value_fp)(shape[0]*shape[1]));
    log_l << "Applying (in-place) FFT on the model image.  Norm = " << norm << LogIO::POST;

    // Below is the following sequence of operations on the GPU: shift -> FFT -> shift.
    // The second shift is equivalent of applying a phase shift corresponding to N/2 pixels. I think...


    hpgGridder_p->shift_model(ShiftDirection::FORWARD);                                 // <============= shift_model()

    err = hpgGridder_p->apply_model_fft(norm,FFTSign::NEGATIVE); // <============= apply_model_fft()
    if (err) log_l << "Gridder::apply_model_fft() failed.  Error type: " << static_cast<int>(err->type()) << LogIO::SEVERE;
    hpgGridder_p->shift_model(ShiftDirection::BACKWARD);                                 // <============= shift_model()

    log_l << "...finished FFT." << LogIO::POST;

    return bool(!err);
  }
  //
  //-----------------------------------------------------------------------------------
  // Template implementation for DataToGrid
  //
  template <class T>
  void AWVisResamplerHPG::DataToGridImpl_p(Array<T>& grid,  VBStore& vbs, 
					Matrix<Double>& sumwt,const Bool& dopsf,
					Bool /*useConjFreqCF*/)
  {
    Int nDataChan, nDataPol, nGridPol, nGridChan, nx, ny, nW;//, nCFFreq;
    Int targetIMPol, rbeg, rend;//, PolnPlane, ConjPlane;
    Int startChan, endChan;
    
    //
    // Load the VB in the vbsList_p.  This list is emptied, once it is
    // full, in the code at the end of this block.
    //

      rbeg = 0;       rend = vbs.nRow_p;
      //cerr << "NROWS " << vbs.nRow_p << endl;
      if(rend==0)
        return;
      rbeg = vbs.beginRow_p;
      rend = vbs.endRow_p;
    
      nx = grid.shape()[0]; ny = grid.shape()[1]; 
      nGridPol = grid.shape()[2]; nGridChan = grid.shape()[3];

      nDataPol  = vbs.flagCube_p.shape()[0];
      nDataChan = vbs.flagCube_p.shape()[1];

      //Bool accumCFs=((vbs.uvw_p.nelements() == 0) && dopsf);
      // Bool Dummy;
      // Double *freq=vbs.freq_p.getStorage(Dummy);
      
      Vector<Double> wVals, fVals; PolMapType mVals, mNdx, conjMVals, conjMNdx;
      Double fIncr, wIncr;
      CountedPtr<CFBuffer> cfb = (*vb2CFBMap_p)[0];
  
     
      // This loads the all-importnat conjMNDx and mNdx maps
      //
      cfb->getCoordList(fVals,wVals,mNdx, mVals, conjMNdx, conjMVals, fIncr, wIncr);
      //    runTimeG1_p += timer_p.real();
      nW = wVals.nelements();

      // timer.mark();
      startChan = 0;
      endChan = nDataChan;
      // if (accumCFs) {startChan = vbs.startChan_p;  endChan = vbs.endChan_p;}
      // else          {startChan = 0;                endChan = nDataChan;}

      Int vbSpw = (vbs.vb_p)->spectralWindows()(0);
      uInt nVisPol=0;// nRows=rend-rbeg+1, nVisChan=endChan - startChan+1;
      for(Int ipol=0; ipol< nDataPol; ipol++)
	{
	  targetIMPol=polMap_p(ipol);
	  if ((targetIMPol>=0) && (targetIMPol<nGridPol)) 
	    nVisPol++;
	}
      //
      // Page-in the new CFs and load them on the device upon SPW change.
      //
      //hpg::CFSimpleIndexer cfsi({1,false},{1,false},{1,true},{1,true}, 1);

      // Reload CFs if the SPW changed and the setup is for AW-Projection (wbAWP=T & nW > 1).
      //bool reloadCFs = (((cachedVBSpw_p != vbSpw) &&                        // when data for a new SPW arrives,
      //			 (vbs.nWPlanes_p > 1)    && (vbs.wbAWP_p==true)) || // if WB A-term and w-term corrections are requested, or
      //			(hpgGridder_p==NULL));                              // if the HPG is uninitialized (first-pass)
      bool reloadCFs=(hpgGridder_p==NULL) || (cachedVBSpw_p != vbSpw);
      //TESTOOO
      //reloadCFs=True;
      //////////////////
      double spwRefFreq = vbs.vb_p->subtableColumns().spectralWindow().refFrequency()(vbSpw);
      int nVBAntenna = vbs.vb_p->nAntennas();
      int nVBChannels = vbs.vb_p->nChannels();
      bool do_degrid=(HPGModelImageName_p != "") || modelImage_p;
      Vector<Int> wNdxList, spwNdxList;
      if (cachedVBSpw_p != vbSpw)
	{
	  cachedVBSpw_p = vbSpw;
	  // LogIO log_l(LogOrigin("AWVisResamplerHPG","DataToGrid_impl"));
	  // log_l << "SPW: " << vbSpw << " " << "Field: " << (vbs.vb_p)->fieldId()(0) << LogIO::POST; 
	  // Set the flag to re-load and re-send the CFs.
	  if (hpgGridder_p==NULL)
	    {
	      // The grid is always a 4D array: NX x NY x NPol x NChan
	      const std::array<unsigned, 4> gridSize{(uInt)nx,(uInt)ny,(uInt)nGridPol,(uInt)nGridChan};
	      //const std::array<float, 2> gridScale{(float)uvwScale_p(0), (float)uvwScale_p(1)};
	      const std::array<double, 2> gridScale{uvwScale_p(0), uvwScale_p(1)};
		      // std::vector<std::array<int, HPGNPOL> > mueller_indexes,conjugate_mueller_indexes;

	      // makeMuellerIndexes<HPGNPOL>(mNdx, mVals, (uInt)nGridPol, mueller_indexes);
	      // makeMuellerIndexes<HPGNPOL>(conjMNdx, conjMVals, (uInt)nGridPol, conjugate_mueller_indexes);
          cerr << "####calling INIT " << " mNdx,mVals,conjMNdx,conjMVals,nVBAntenna, nVBChannels " <<  mNdx << " " << mVals << " " << conjMNdx << " " << conjMVals << " " << nVBAntenna << " " << nVBChannels << endl;
	      hpgGridder_p = initGridder2<HPGNPOL>(HPGDevice_p,1,&cfArray, gridSize, gridScale,mNdx,mVals,conjMNdx,conjMVals,
						   nVBAntenna, nVBChannels);
						   //mueller_indexes,conjugate_mueller_indexes);
	      do_degrid = sendModelImage(gridSize);
	      if (do_degrid)
		{
		  LogIO log_l(LogOrigin("AWVisResamplerHPG","DataToGrid_impl"));
		  log_l << "Running degrid_grid GPU kernel" << LogIO::POST;
		}
		
	    }

	  if (reloadCFs)
	    {
	      // If wplanes == 1, make a list of all SPW CFs.  If
	      // wplanes > 1, make a list of all W CFs for the current
	      // SPW.
	      //int tspw=(vbs.nWPlanes_p == 1)?-1:vbSpw;

	      makeAWLists(*cfb, vbs.wbAWP_p, vbs.nWPlanes_p,
			  vbs.imRefFreq_p, spwRefFreq,
			  wNdxList, spwNdxList,
                          vbSpw);

	      // MakeCFArray::makeCFArray(vbs.wbAWP_p, vbs.nWPLanes_p,
	      // 		  ...(skyImage),...
	      //               polMap_p,
	      // 		  ...(cfs2_l),
	      //               tswp,//vbSpw_l
	      // 		  spwRefFreq, nDataPol, nGridPol, wNdxList,spwNdxList)
	      // MakeCFArray::makeCFArray_p(cfArray,
	      // 				 vbs.paQuant_p.getValue("rad"),// vbPA,
	      // 				 vbs.imRefFreq_p,
	      // 				 *cfb,
	      // 				 nGridPol,
	      // 				 nDataPol,
	      // 				 polMap_p,
	      // 				 wNdxList,
	      // 				 spwNdxList);


	      //cerr << "SPW ID: " << (vbs.vb_p)->spectralWindows()[0] << endl;
	      //cerr << "ImRefFreq, SpwRefFreq: " << vbs.imRefFreq_p << " " << spwRefFreq << endl;
	      //cerr << "WList: "; for(auto& n : wNdxList) cerr << n << ' '; cerr << endl;
	      //cerr << "SPWList: "; for(auto& n : spwNdxList) cerr << n << ' '; cerr << endl;

	      //cfsi = loadCF(vbs, nGridPol, nDataPol,true);

	      // Returns std::tuple<hpg::opt_t<hpg::Error>, hpg::CFSimpleIndexer>
	      auto ret= loadCF(vbs, nGridPol, nDataPol,wNdxList,spwNdxList,true);
	      //auto ret= loadCF(cfb,vbs, nGridPol, nDataPol,wNdxList,spwNdxList,true);
	      cfsi_p=std::get<1>(ret);
	      if (std::get<0>(ret))
		{
		  LogIO log_l(LogOrigin("AWVisResamplerHPG","DataToGrid_impl"));
		  log_l << "HPGError in calling set_convolution_function(): "  
			<< " " << static_cast<int>(std::get<0>(ret)->type()) << LogIO::EXCEPTION;
		}
	    }
	}
      
      //cerr << "FIELD ID: " << (vbs.vb_p)->fieldId()[0] << " ";


     
      casacore::Vector<casacore::Vector<casacore::Double> > pointingOffsets= cfb->getPointingOffset();
      //cerr << "SHAPE of PTING off "<<  pointingOffsets.shape() <<  " max " << pointingOffsets[0] << endl;
      // rbeg=19;
      // rend=20; //For a single-vis test
      //startChan=32; endChan=33;
      //      cerr << chanMap_p << endl;

      //cerr << "Pointing offsets: " << pointingOffsets << " " << rbeg << " " << rend << " " << startChan << " " << endChan << " " << max(vbs.vb_p->visCube()) << endl;
      std::vector<hpg::VisData<HPGNPOL> >
	hpgVB =                                                       makeHPGVisBuffer<HPGNPOL>(sumwt,
					  vbs,
					  vb2CFBMap_p,
					  nGridPol, nGridChan, nVisPol,
					  rbeg, rend, startChan, endChan,
					  chanMap_p, polMap_p,
					  conjMNdx, conjMVals,
					  dphase_p,pointingOffsets,
					  cfsi_p);

      // Add to the list only if the hpgVB holds any data. This guard
      // is required since in the loop below we also count the number
      // visibilities gridded (the nVisGridded_p).
      if (hpgVB.size() > 0) hpgVBList_p.push_back(hpgVB);
      //
      // If the vbsList_p is full, send the VBs loaded in vbsList_p for gridding, and empty vbsList_p.
      // std::move() empties the storage of its argument.
      //
      timer_p.mark();
      cerr << "hpgVBList_p.size  " <<  hpgVBList_p.size() << " maxVBList " <<  maxVBList_p << endl;
      if (hpgVBList_p.size() >= maxVBList_p)
	{
	  for(unsigned i=0;i<hpgVBList_p.size();i++)
	    {
	      nVisGridded_p += hpgVBList_p[i].size()*HPGNPOL;
	      hpg::opt_t<hpg::Error> err;
	      if (do_degrid)
		err = hpgGridder_p->degrid_grid_visibilities(
							     hpg::Device::OpenMP,
							     std::move(hpgVBList_p[i])
							     );
	      else
		err = hpgGridder_p->grid_visibilities(
						      hpg::Device::OpenMP,
						      std::move(hpgVBList_p[i])
						      );

	      if (err)
		{
		  LogIO log_l(LogOrigin("AWVisResamplerHPG","DataToGrid_impl"));
		  log_l << "Failed hpg::" << ((do_degrid)?"degrid_grid_visibilities()":"grid_visibilities()") << " Error type: " << static_cast<int>(err->type()) << LogIO::SEVERE;
		}
	    }
	  hpgVBList_p.resize(0);
	}
      griddingTime += timer_p.real();
      return;
  }
  //
  //-----------------------------------------------------------------------------------
  //
  template <class T>
  void AWVisResamplerHPG::DataToGridImpl2_p(Array<T>& grid,  VBStore& vbs,
                    AWConvFuncHolder& awh, 
					Matrix<Double>& sumwt,const Bool& dopsf)
  {
    Int nDataChan, nDataPol, nGridPol, nGridChan, nx, ny, nW;//, nCFFreq;
    Int targetIMPol, rbeg, rend;//, PolnPlane, ConjPlane;
    Int startChan, endChan;
    

      rbeg = 0;       rend = vbs.nRow_p;
      //cerr << "NROWS " << vbs.nRow_p << endl;
      if(rend==0)
        return;
      rbeg = vbs.beginRow_p;
      rend = vbs.endRow_p;
    
      nx = grid.shape()[0]; ny = grid.shape()[1]; 
      nGridPol = grid.shape()[2]; nGridChan = grid.shape()[3];

      nDataPol  = vbs.flagCube_p.shape()[0];
      nDataChan = vbs.flagCube_p.shape()[1];

      // timer.mark();
      startChan = 0;
      endChan = nDataChan;

      Int vbSpw = (vbs.vb_p)->spectralWindows()(0);
      uInt nVisPol=0;// nRows=rend-rbeg+1, nVisChan=endChan - startChan+1;
      for(Int ipol=0; ipol< nDataPol; ipol++)
	{
	  targetIMPol=polMap_p(ipol);
	  if ((targetIMPol>=0) && (targetIMPol<nGridPol)) 
	    nVisPol++;
	}
      //
   
      bool reloadCFs=(hpgGridder_p==NULL) || (cachedVBSpw_p != vbSpw);
      //bool reloadCFs=(hpgGridder_p==NULL);
      //TESTOOO
      //reloadCFs=True;
      //////////////////
      //cerr <<  "vbspw" <<  vbSpw <<  "  cached " <<  cachedVBSpw_p <<  endl;
      
      double spwRefFreq = vbs.vb_p->subtableColumns().spectralWindow().refFrequency()(vbSpw);
      int nVBAntenna = vbs.vb_p->nAntennas();
      int nVBChannels = vbs.vb_p->nChannels();
      bool do_degrid=(HPGModelImageName_p != "") || modelImage_p;
      Vector<Int> wNdxList, spwNdxList;
      if (cachedVBSpw_p != vbSpw)
	{
	  cachedVBSpw_p = vbSpw;
	  // LogIO log_l(LogOrigin("AWVisResamplerHPG","DataToGrid_impl"));
	  // log_l << "SPW: " << vbSpw << " " << "Field: " << (vbs.vb_p)->fieldId()(0) << LogIO::POST; 
	  // Set the flag to re-load and re-send the CFs.
	  if (hpgGridder_p==NULL)
	    {
	      // The grid is always a 4D array: NX x NY x NPol x NChan
	      const std::array<unsigned, 4> gridSize{(uInt)nx,(uInt)ny,(uInt)nGridPol,(uInt)nGridChan};
	      //const std::array<float, 2> gridScale{(float)uvwScale_p(0), (float)uvwScale_p(1)};
	      const std::array<double, 2> gridScale{uvwScale_p(0), uvwScale_p(1)};
		      // std::vector<std::array<int, HPGNPOL> > mueller_indexes,conjugate_mueller_indexes;

	      // makeMuellerIndexes<HPGNPOL>(mNdx, mVals, (uInt)nGridPol, mueller_indexes);
	      // makeMuellerIndexes<HPGNPOL>(conjMNdx, conjMVals, (uInt)nGridPol, conjugate_mueller_indexes);
          //cerr << "####calling INIT3 " << " nVBAntenna, nVBChannels " << " " << nVBAntenna << " " << nVBChannels << endl;
	      hpgGridder_p = initGridder3<HPGNPOL>(HPGDevice_p,1,&cfArray, gridSize, gridScale,nVBAntenna, nVBChannels);
						   //mueller_indexes,conjugate_mueller_indexes);
	      do_degrid = sendModelImage(gridSize);
	      if (do_degrid)
		{
		  LogIO log_l(LogOrigin("AWVisResamplerHPG","DataToGrid_impl"));
		  log_l << "Running degrid_grid GPU kernel" << LogIO::POST;
		}
		
	    }

	  if (reloadCFs)
	    {
	      



	      
	      auto ret= loadCF(vbs, awh, true);
	      //auto ret= loadCF(cfb,vbs, nGridPol, nDataPol,wNdxList,spwNdxList,true);
	      cfsi_p=std::get<1>(ret);
	      if (std::get<0>(ret))
		{
		  LogIO log_l(LogOrigin("AWVisResamplerHPG","DataToGrid_impl"));
		  log_l << "HPGError in calling set_convolution_function(): "  
			<< " " << static_cast<int>(std::get<0>(ret)->type()) << LogIO::EXCEPTION;
		}
	    }
	}
      
      //cerr << "FIELD ID: " << (vbs.vb_p)->fieldId()[0] << " ";


     
      //casacore::Vector<casacore::Vector<casacore::Double> > pointingOffsets(1);
      //cerr << "SHAPE of PTING off "<<  pointingOffsets.shape() <<  " max " << pointingOffsets[0] << endl;
      // rbeg=19;
      // rend=20; //For a single-vis test
      //startChan=32; endChan=33;
      //      cerr << chanMap_p << endl;

      //cerr << "Pointing offsets: " << pointingOffsets << " " << rbeg << " " << rend << " " << startChan << " " << endChan << " " << max(vbs.vb_p->visCube()) << endl;
      std::vector<hpg::VisData<HPGNPOL> > hpgVB = 
		makeHPGVisBuffer2<HPGNPOL>(sumwt, vbs,awh,nGridPol, nGridChan, rbeg, rend,startChan, endChan,chanMap_p,polMap_p, dphase_p,cfsi_p);
      // Add to the list only if the hpgVB holds any data. This guard
      // is required since in the loop below we also count the number
      // visibilities gridded (the nVisGridded_p).
      if (hpgVB.size() > 0) hpgVBList_p.push_back(hpgVB);
      //
      // If the vbsList_p is full, send the VBs loaded in vbsList_p for gridding, and empty vbsList_p.
      // std::move() empties the storage of its argument.
      //
      timer_p.mark();
      //cerr << "hpgVBList_p.size  " <<  hpgVBList_p.size() << " maxVBList " <<  maxVBList_p << endl;
      if (hpgVBList_p.size() >= maxVBList_p)
	{
	  for(unsigned i=0;i<hpgVBList_p.size();i++)
	    {
	      nVisGridded_p += hpgVBList_p[i].size()*HPGNPOL;
	      hpg::opt_t<hpg::Error> err;
	      if (do_degrid)
		err = hpgGridder_p->degrid_grid_visibilities(
							     hpg::Device::OpenMP,
							     std::move(hpgVBList_p[i])
							     );
	      else
		err = hpgGridder_p->grid_visibilities(
						      hpg::Device::OpenMP,
						      std::move(hpgVBList_p[i])
						      );

	      if (err)
		{
		  LogIO log_l(LogOrigin("AWVisResamplerHPG","DataToGrid_impl"));
		  log_l << "Failed hpg::" << ((do_degrid)?"degrid_grid_visibilities()":"grid_visibilities()") << " Error type: " << static_cast<int>(err->type()) << LogIO::SEVERE;
		}
	    }
	  hpgVBList_p.resize(0);
	}
      griddingTime += timer_p.real();
      return;
  }
  //
  template
  void AWVisResamplerHPG::DataToGridImpl2_p(Array<DComplex>& grid, VBStore& vbs, AWConvFuncHolder& awh, 
  					Matrix<Double>& sumwt,const Bool& dopsf); 
  template
  void AWVisResamplerHPG::DataToGridImpl2_p(Array<Complex>& grid, VBStore& vbs, AWConvFuncHolder& awh, 
					Matrix<Double>& sumwt,const Bool& dopsf); 
  //-----------------------------------------------------------------------------------
  //
  std::shared_ptr<std::complex<double>> AWVisResamplerHPG::getGridPtr(size_t& size) const
  {
    if (hpgGridder_p)
      {
        hpgGridder_p->fence();
        size = hpgGridder_p->grid_values_span();
        return hpgGridder_p->grid_values_ptr();
      }
    else
      {
        size = 0;
        return std::shared_ptr<std::complex<double>>();
      }
  }

  std::shared_ptr<double> AWVisResamplerHPG::getSumWeightsPtr(size_t& size) const
  {
    if (hpgGridder_p)
      {
        hpgGridder_p->fence();
        size = hpgGridder_p->grid_weights_span();
        return hpgGridder_p->grid_weights_ptr();
      }
    else
      {
        size = 0;
        return std::shared_ptr<double>();
      }
  }


  using namespace casacore;
};// end namespace casa
#endif  // endif of USE_HPG
