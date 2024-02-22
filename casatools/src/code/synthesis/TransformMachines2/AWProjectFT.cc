// -*- C++ -*-
//# AWProjectFT.cc: Implementation of AWProjectFT class
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

#include <casacore/casa/Quanta/UnitMap.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/casa/Quanta/UnitVal.h>
#include <casacore/casa/Containers/Block.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/OS/HostInfo.h>
#include <sstream>

#include <casacore/coordinates/Coordinates/CoordinateSystem.h>
#include <casacore/images/Images/ImageInterface.h>

#include <synthesis/TransformMachines/StokesImageUtil.h>
#include <synthesis/TransformMachines/SynthesisError.h>
#include <synthesis/TransformMachines2/AWProjectFT.h>
#include <synthesis/TransformMachines2/CFStore2.h>
#include <synthesis/MeasurementComponents/ExpCache.h>
#include <synthesis/MeasurementComponents/CExp.h>
#include <synthesis/TransformMachines2/AWVisResampler.h>
#include <synthesis/TransformMachines2/VBStore.h>

#include <casacore/scimath/Mathematics/FFTServer.h>
#include <casacore/scimath/Mathematics/MathFunc.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <iostream>
#include <casacore/casa/OS/Timer.h>

#include <synthesis/TransformMachines2/ATerm.h>
#include <synthesis/TransformMachines2/NoOpATerm.h>
#include <synthesis/TransformMachines2/AWConvFunc.h>
#include <synthesis/TransformMachines2/EVLAAperture.h>
#include <iomanip>
//#define CONVSIZE (1024*2)
// #define OVERSAMPLING 2
#define USETABLES 0           // If equal to 1, use tabulated exp() and
			      // complex exp() functions.
#define MAXPOINTINGERROR 250.0 // Max. pointing error in arcsec used to
// determine the resolution of the
// tabulated exp() function.
#define DORES true


using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN
  using namespace vi;  
#define NEED_UNDERSCORES
  namespace refim{
  extern "C" 
  {
    //
    // The Gridding Convolution Function (GCF) used by the underlying
    // gridder written in FORTRAN.
    //
    // The arguments must all be pointers and the value of the GCF at
    // the given (u,v) point is returned in the weight variable.  Making
    // this a function which returns a complex value (namely the weight)
    // has problems when called in FORTRAN - I (SB) don't understand
    // why.
    //
#if defined(NEED_UNDERSCORES)
#define nwcppeij nwcppeij_
#endif
    //
    //---------------------------------------------------------------
    //
    IlluminationConvFunc awEij2;
    void awcppeij2(Double *griduvw, Double *area,
		 Double *raoff1, Double *decoff1,
		 Double *raoff2, Double *decoff2, 
		 Int *doGrad,
		 Complex *weight,
		 Complex *dweight1,
		 Complex *dweight2,
		 Double *currentCFPA)
    {
      Complex w,d1,d2;
      awEij2.getValue(griduvw, raoff1, raoff2, decoff1, decoff2,
		    area,doGrad,w,d1,d2,*currentCFPA);
      *weight   = w;
      *dweight1 = d1;
      *dweight2 = d2;
    }
  }

    ATerm* AWProjectFT::createTelescopeATerm(const String& telescopeName,
					     const Bool& // isATermOn
					     )
  {
    
    //    if (!isATermOn) return new NoOpATerm();
    
    // MSObservationColumns msoc(ms.observation());
    // String ObsName=msoc.telescopeName()(0);
    if ((telescopeName == "EVLA") || (telescopeName == "VLA"))
      return new EVLAAperture();
    else
      {
	LogIO os(LogOrigin("AWProjectFT2", "createTelescopeATerm",WHERE));
	os << "Telescope name ('"+
	  telescopeName+"') in the MS not recognized to create the telescope specific ATerm" 
	   << LogIO::EXCEPTION;
      }
    
    return NULL;
  }

  CountedPtr<ConvolutionFunction> AWProjectFT::makeCFObject(const String& telescopeName,
							    const Bool aTermOn,
							    const Bool psTermOn,
							    const Bool wTermOn,
							    const Bool,// mTermOn,
							    const Bool wbAWP,
							    const Bool conjBeams)
  {
    (void)wTermOn;//Unused
    CountedPtr<ATerm> apertureFunction = AWProjectFT::createTelescopeATerm(telescopeName, aTermOn);
    CountedPtr<PSTerm> psTerm = new PSTerm();
    CountedPtr<WTerm> wTerm = new WTerm();
    //if (cfBufferSize > 0) apertureFunction->setConvSize(cfBufferSize);
    //
    // Selectively switch off CFTerms.
    //
    if (aTermOn == false) {apertureFunction->setOpCode(CFTerms::NOOP);}
    if (psTermOn == false) psTerm->setOpCode(CFTerms::NOOP);
    if (wTermOn == False) wTerm->setOpCode(CFTerms::NOOP);
    //
    // Construct the CF object with appropriate CFTerms.
    //
    CountedPtr<ConvolutionFunction> awConvFunc;
    //    awConvFunc = new AWConvFunc(apertureFunction,psTerm,wTerm, !wbAWP);
    //if ((ftmName=="mawprojectft") || (mTermOn))

    awConvFunc = new AWConvFunc(apertureFunction,psTerm,wTerm,wbAWP, conjBeams);

    return awConvFunc;
  }
  //
  //---------------------------------------------------------------
  //
  AWProjectFT::AWProjectFT()
    : FTMachine(), padding_p(1.0), nWPlanes_p(1),
      imageCache(0), cachesize(0), tilesize(16),
      gridder(0), isTiled(false),  lattice( ), 
      maxAbsData(0.0), centerLoc(IPosition(4,0)), offsetLoc(IPosition(4,0)),
      pointingToImage(0), usezero_p(false), avgPB_p(nullptr), 
      epJ_p(nullptr),
      doPBCorrection(true), conjBeams_p(true),/*cfCache_p(cfcache),*/ paChangeDetector(),
      rotateOTFPAIncr_p(0.1),
      Second("s"),Radian("rad"),Day("d"), pbNormalized_p(false), paNdxProcessed_p(),
      visResampler_p(nullptr), sensitivityPatternQualifier_p(-1),sensitivityPatternQualifierStr_p(""),
      rotatedConvFunc_p(),
      runTime1_p(0.0), previousSPWID_p(-1), self_p(nullptr), vb2CFBMap_p(nullptr), po_p(nullptr),wbAWP_p(true),
    timemass_p(0.0), timegrid_p(0.0), timedegrid_p(0.0)
  {
    //    convSize=0;
    tangentSpecified_p=false;
    lastIndex_p=0;
    paChangeDetector.reset();
    pbLimit_p=5e-2;
    //
    // Get various parameters from the visibilities.  
    //
    doPointing=1; 

    maxConvSupport=-1;  
    //
    // Set up the Conv. Func. disk cache manager object.
    //
    // if (!cfCache_p.null()) delete &cfCache_p;
    // cfCache_p=cfcache;
    convSampling=-1;
    //convSize=CONVSIZE;
    Long hostRAM = (HostInfo::memoryTotal(true)*1024); // In bytes
    hostRAM = hostRAM/(sizeof(Float)*2); // In complex pixels
    if (cachesize > hostRAM) cachesize=hostRAM;
    sigma=1.0;
    canComputeResiduals_p=DORES;
    // cfs2_p = &cfCache_p->memCache2_p[0];//new CFStore2;
    // cfwts2_p =  &cfCache_p->memCacheWt2_p[0];//new CFStore2;
    pop_p->init();
    CFBuffer::initCFBStruct(cfbst_pub);
    //    rotatedConvFunc_p.data=new Array<Complex>();    
    //    self_p.reset(this);
    vb2CFBMap_p = new VB2CFBMap();
    po_p = new PointingOffsets();
    po_p->setOverSampling(convSampling);
  }
  //
  //---------------------------------------------------------------
  //
  AWProjectFT::AWProjectFT(Int nWPlanes, Long icachesize, 
			   CountedPtr<CFCache>& cfcache,
			   CountedPtr<ConvolutionFunction>& cf,
			   CountedPtr<VisibilityResamplerBase>& visResampler,
			   Bool applyPointingOffset,
			   vector<float> pointingOffsetSigDev,
			   Bool doPBCorr,
			   Int itilesize, 
			   Float pbLimit,
			   Bool usezero,
			   Bool conjBeams,
			   Bool doublePrecGrid,
			   PolOuterProduct::MuellerType muellerType)
    : FTMachine(cfcache,cf), padding_p(1.0), nWPlanes_p(nWPlanes),
      imageCache(0), cachesize(icachesize), tilesize(itilesize),
      gridder(0), isTiled(false),  lattice( ), 
      maxAbsData(0.0), centerLoc(IPosition(4,0)), offsetLoc(IPosition(4,0)),
      pointingToImage(0), usezero_p(usezero), avgPB_p(nullptr),
      // convFunc_p(), convWeights_p(),
      epJ_p(),
      doPBCorrection(doPBCorr), conjBeams_p(conjBeams), 
      /*cfCache_p(cfcache),*/ paChangeDetector(),
      rotateOTFPAIncr_p(0.1),
      Second("s"),Radian("rad"),Day("d"), pbNormalized_p(false),
      visResampler_p(visResampler), sensitivityPatternQualifier_p(-1),sensitivityPatternQualifierStr_p(""),
      rotatedConvFunc_p(), runTime1_p(0.0),  previousSPWID_p(-1),self_p(nullptr), vb2CFBMap_p(nullptr), po_p(nullptr),wbAWP_p(true), timemass_p(0.0), timegrid_p(0.0), timedegrid_p(0.0)
  {
    //convSize=0;
    tangentSpecified_p=false;
    lastIndex_p=0;
    paChangeDetector.reset();
    pbLimit_p=pbLimit;
    //
    // Get various parameters from the visibilities.  
    //
    if (applyPointingOffset) doPointing=1; else doPointing=0;

    maxConvSupport=-1;  
    //
    // Set up the Conv. Func. disk cache manager object.
    //
    // if (!cfCache_p.null()) delete &cfCache_p;
    // cfCache_p=cfcache;
    convSampling=convFuncCtor_p->getOversampling();
    //convSize=CONVSIZE;
    Long hostRAM = (HostInfo::memoryTotal(true)*1024); // In bytes
    hostRAM = hostRAM/(sizeof(Float)*2); // In complex pixels
    if (cachesize > hostRAM) cachesize=hostRAM;
    sigma=1.0;
    canComputeResiduals_p=DORES;
    if (!cfCache_p.null())
      {
	cfs2_p = CountedPtr<CFStore2>(&(cfCache_p->memCache2_p)[0],false);//new CFStore2;
	cfwts2_p =  CountedPtr<CFStore2>(&cfCache_p->memCacheWt2_p[0],false);//new CFStore2;
      }
    pop_p->init();
    useDoubleGrid_p=doublePrecGrid;
    //    rotatedConvFunc_p.data=new Array<Complex>();
    CFBuffer::initCFBStruct(cfbst_pub);
    muellerType_p = muellerType;
    //    self_p.reset(this);
    vb2CFBMap_p = new VB2CFBMap();
    po_p = new PointingOffsets();
    po_p->setOverSampling(convSampling);
    vb2CFBMap_p->setPOSigmaDev(pointingOffsetSigDev);
    wbAWP_p=convFuncCtor_p->isWBAWP();

  }
  //
  //---------------------------------------------------------------
  //
  AWProjectFT::AWProjectFT(const RecordInterface& stateRec)
    : FTMachine(),Second("s"),Radian("rad"),Day("d"),visResampler_p(nullptr), self_p(nullptr), vb2CFBMap_p(nullptr), po_p(nullptr),wbAWP_p(true), timemass_p(0.0), timegrid_p(0.0), timedegrid_p(0.0)
  {
    //
    // Construct from the input state record
    //
    String error;
    
    if (!fromRecord(stateRec)) {
      LogIO log_l(LogOrigin("AWProjectFT2", "AWProjectFT[R&D]"));
      log_l << "Failed to create " << name() << " object." << LogIO::EXCEPTION;
    };
    maxConvSupport=-1;
    convSampling=-1;
    visResampler_p->init(useDoubleGrid_p);
    //convSize=CONVSIZE;
    canComputeResiduals_p=DORES;
    if (!cfCache_p.null())
      {
	cfs2_p = CountedPtr<CFStore2>(&cfCache_p->memCache2_p[0],false);//new CFStore2;
	cfwts2_p =  CountedPtr<CFStore2>(&cfCache_p->memCacheWt2_p[0],false);//new CFStore2;
      }
    pop_p->init();
    //    self_p.reset(this);
    vb2CFBMap_p = new VB2CFBMap();
    po_p = new PointingOffsets();
    po_p->setOverSampling(convSampling);
  }
  //
  //----------------------------------------------------------------------
  //
  AWProjectFT::AWProjectFT(const AWProjectFT& other):FTMachine()
  {
    operator=(other);
  }
  //
  //---------------------------------------------------------------
  //
  // This is nasty, we should use CountedPointers here.
  AWProjectFT::~AWProjectFT() 
  {
      if(imageCache) delete imageCache;
      imageCache=0;
      if(gridder) delete gridder;
      gridder=0;
  }
  //
  //---------------------------------------------------------------
  //
  AWProjectFT& AWProjectFT::operator=(const AWProjectFT& other)
  {
    if(this!=&other) 
      {
	//Do the base parameters
	FTMachine::operator=(other);

	
	padding_p=other.padding_p;
	
	nWPlanes_p=other.nWPlanes_p;
	imageCache=other.imageCache;
	cachesize=other.cachesize;
	tilesize=other.tilesize;
	cfRefFreq_p = other.cfRefFreq_p;
	if(other.gridder==0)
	  gridder=0;
	else
	  {
	    uvScale.resize();
	    uvOffset.resize();
	    uvScale=other.uvScale;
	    uvOffset=other.uvOffset;
	    gridder = new ConvolveGridder<Double, Complex>(IPosition(2, nx, ny),
							   uvScale, uvOffset,
							   "SF");
	  }

	isTiled=other.isTiled;
	lattice=0;

	maxAbsData=other.maxAbsData;
	centerLoc=other.centerLoc;
	offsetLoc=other.offsetLoc;
	pointingToImage=other.pointingToImage;
	usezero_p=other.usezero_p;

	
	padding_p=other.padding_p;
	imageCache=other.imageCache;
	cachesize=other.cachesize;
	tilesize=other.tilesize;
	isTiled=other.isTiled;
	maxAbsData=other.maxAbsData;
	centerLoc=other.centerLoc;
	offsetLoc=other.offsetLoc;
	pointingToImage=other.pointingToImage;
	 useDoubleGrid_p=other.useDoubleGrid_p;
	doPBCorrection = other.doPBCorrection;
	maxConvSupport= other.maxConvSupport;

	epJ_p=other.epJ_p;
	//convSize=other.convSize;
	lastIndex_p=other.lastIndex_p;
	paChangeDetector=other.paChangeDetector;
	pbLimit_p=other.pbLimit_p;
	//
	// Get various parameters from the visibilities.  
	//
	doPointing=other.doPointing;

	maxConvSupport=other.maxConvSupport;
	//
	// Set up the Conv. Func. disk cache manager object.
	//
	cfCache_p=other.cfCache_p;
	convSampling=other.convSampling;
	//convSize=other.convSize;
	cachesize=other.cachesize;
    
	currentCFPA=other.currentCFPA;
	lastPAUsedForWtImg = other.lastPAUsedForWtImg;
	avgPB_p = other.avgPB_p;

	convFuncCtor_p = other.convFuncCtor_p;
	pbNormalized_p = other.pbNormalized_p;
	sensitivityPatternQualifier_p = other.sensitivityPatternQualifier_p;
	sensitivityPatternQualifierStr_p = other.sensitivityPatternQualifierStr_p;
	visResampler_p=other.visResampler_p; // Copy the counted pointer
	//	visResampler_p=other.visResampler_p->clone();
	//	*visResampler_p = *other.visResampler_p; // Call the appropriate operator=()

	rotatedConvFunc_p = other.rotatedConvFunc_p;
	cfs2_p = other.cfs2_p;
	cfwts2_p = other.cfwts2_p;
	paNdxProcessed_p = other.paNdxProcessed_p;
	imRefFreq_p = other.imRefFreq_p;
	conjBeams_p = other.conjBeams_p;
	rotateOTFPAIncr_p=other.rotateOTFPAIncr_p;
	computePAIncr_p=other.computePAIncr_p;
	runTime1_p = other.runTime1_p;
	muellerType_p = other.muellerType_p;
	previousSPWID_p = other.previousSPWID_p;
	vb2CFBMap_p = other.vb2CFBMap_p;
	po_p = other.po_p;
	//	self_p = other.self_p;
	wbAWP_p=other.wbAWP_p;
        timemass_p=0.0;
        timegrid_p = 0.0;
        timedegrid_p = 0.0;
      };
    return *this;
  };
  //
  //----------------------------------------------------------------------
  //
  void AWProjectFT::init(const vi::VisBuffer2& /*vb*/) 
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
    
    uvOffset.resize(3);
    uvOffset(0)=nx/2;
    uvOffset(1)=ny/2;
    uvOffset(2)=0;
    
    if(gridder) delete gridder;
    gridder=0;
    gridder = new ConvolveGridder<Double, Complex>(IPosition(2, nx, ny),
						   uvScale, uvOffset,
						   "SF");
    
    // Set up image cache needed for gridding. 
    if(imageCache) delete imageCache;
    imageCache=0;
    
    // The tile size should be large enough that the
    // extended convolution function can fit easily
    if(isTiled) 
      {
	Float tileOverlap=0.5;
	tilesize=min(256,tilesize);
	IPosition tileShape=IPosition(4,tilesize,tilesize,npol,nchan);
	Vector<Float> tileOverlapVec(4);
	tileOverlapVec=0.0;
	tileOverlapVec(0)=tileOverlap;
	tileOverlapVec(1)=tileOverlap;
	if (sizeof(long) < 4)  // 32-bit machine
	  {
	    Int tmpCacheVal=static_cast<Int>(cachesize);
	    imageCache=new LatticeCache <Complex> (*image, tmpCacheVal, tileShape, 
						   tileOverlapVec,
						   (tileOverlap>0.0));
	  }
	else  // 64-bit machine
	  {
	    Long tmpCacheVal=cachesize;
	    imageCache=new LatticeCache <Complex> (*image, tmpCacheVal, tileShape, 
						   tileOverlapVec,
						   (tileOverlap>0.0));
	  }
      }
    
#if(USETABLES)
    Double StepSize;
    Int N=500000;
    StepSize = abs((((2*nx)/uvScale(0))/(sigma) + 
		    MAXPOINTINGERROR*1.745329E-02*(sigma)/3600.0))/N;
    if (!awEij2.isReady())
      {
	log_l << "Making lookup table for exp function with a resolution of " 
	      << StepSize << " radians.  "
	      << "Memory used: " << sizeof(Float)*N/(1024.0*1024.0)<< " MB." 
	      << LogIO::NORMAL 
	      <<LogIO::POST;
	
	awEij2.setSigma(sigma);
	awEij2.initExpTable(N,StepSize);
	//    ExpTab.build(N,StepSize);
	
	log_l << "Making lookup table for complex exp function with a resolution of " 
	      << 2*M_PI/N << " radians.  "
	      << "Memory used: " << 2*sizeof(Float)*N/(1024.0*1024.0) << " MB." 
	      << LogIO::NORMAL
	      << LogIO::POST;
	awEij2.initCExpTable(N);
	//    CExpTab.build(N);
      }
#endif
    //    vpSJ->reset();
    paChangeDetector.reset();
    makingPSF = false;
    
    log_l << "Using " << visResampler_p->name() << " visibility resampler (a.k.a. \"gridder/degridder\")" << LogIO::POST;
    vb2CFBMap_p->computePhaseScreen_p = visResampler_p->needCFPhaseScreen();
    //cerr << "Current runTime = " << runTime << endl;
  }
  //
  //---------------------------------------------------------------
  //
  MDirection::Convert AWProjectFT::makeCoordinateMachine(const VisBuffer2& vb,
							 const MDirection::Types& From,
							 const MDirection::Types& To,
							 MEpoch& last)
  {
    LogIO log_l(LogOrigin("AWProjectFT2","makeCoordinateMachine[R&D]"));
    Double time = getCurrentTimeStamp(vb);
    
    MEpoch epoch(Quantity(time,Second),MEpoch::TAI);
    //  epoch = MEpoch(Quantity(time,Second),MEpoch::TAI);
    //
    // ...now make an object to hold the observatory position info...
    //
    MPosition pos;
    String ObsName=(vb.subtableColumns()).observation().telescopeName()(vb.arrayId()(0));
    
    if (!MeasTable::Observatory(pos,ObsName))
      log_l << "Observatory position for "+ ObsName + " not found"
	    << LogIO::EXCEPTION;
    //
    // ...now make a Frame object out of the observatory position and
    // time objects...
    //
    MeasFrame frame(epoch,pos);
    //
    // ...finally make the convert machine.
    //
    MDirection::Convert mac(MDirection::Ref(From,frame),
			    MDirection::Ref(To,frame));
    
    MEpoch::Convert toLAST = MEpoch::Convert(MEpoch::Ref(MEpoch::TAI,frame),
					     MEpoch::Ref(MEpoch::LAST,frame));
    last = toLAST(epoch);
    
    return mac;
  }
  //
  //---------------------------------------------------------------
  //
  int AWProjectFT::findPointingOffsets(const VisBuffer2& vb, 
					Array<Float> &l_off,
					Array<Float> &m_off,
					Bool Evaluate)
  {
    //    LogIO log_l(LogOrigin("AWProjectFT2", "findPointingOffsets[R&D]"));
    Int NAnt = 0;
    MEpoch LAST;
    Double thisTime = getCurrentTimeStamp(vb);
    //    Array<Float> pointingOffsets = epJ->nearest(thisTime);
    if (epJ_p.null()) return 0;
    Array<Float> pointingOffsets; epJ_p->nearest(thisTime,pointingOffsets);
    NAnt=pointingOffsets.shape()(2);
    l_off.resize(IPosition(3,1,1,NAnt)); // Poln x NChan x NAnt 
    m_off.resize(IPosition(3,1,1,NAnt)); // Poln x NChan x NAnt 
    // Can't figure out how to do the damn slicing of [Pol,NChan,NAnt,1] array
    // into [Pol,NChan,NAnt] array
    //
    IPosition tndx(3,0,0,0), sndx(4,0,0,0,0);
    for(tndx(2)=0;tndx(2)<NAnt; tndx(2)++,sndx(2)++)
      {
	sndx(0)=0; l_off(tndx) = pointingOffsets(sndx);
	sndx(0)=2; m_off(tndx) = pointingOffsets(sndx);
      }
    return NAnt;
    if (!Evaluate) return NAnt;
    
    //
    // Make a Coordinate Conversion Machine to go from (Az,El) to
    // (HA,Dec).
    //
    MDirection::Convert toAzEl = makeCoordinateMachine(vb,MDirection::HADEC,
						       MDirection::AZEL,
						       LAST);
    MDirection::Convert toHADec = makeCoordinateMachine(vb,MDirection::AZEL,
							MDirection::HADEC,
							LAST);
    //
    // ...and now hope that it all works and works correctly!!!
    //
    Quantity dAz(0,Radian),dEl(0,Radian);
    //
    // An array of shape [2,1,1]!
    //
    Array<Double> phaseDir = vb.subtableColumns().field().phaseDir().getColumn();
    Double RA0   = phaseDir(IPosition(3,0,0,0));
    Double Dec0  = phaseDir(IPosition(3,1,0,0));
    //  
    // Compute reference (HA,Dec)
    //
    Double LST   = LAST.get(Day).getValue();
    Double SDec0 = sin(Dec0), CDec0=cos(Dec0);
    LST -= floor(LST); // Extract the fractional day
    LST *= 2*C::pi;// Convert to Raidan
    
    Double HA0;
    HA0 = LST - RA0;
    Quantity QHA0(HA0,Radian), QDEC0(Dec0,Radian);
    //
    // Convert reference (HA,Dec) to reference (Az,El)
    //
    MDirection PhaseCenter(QHA0, QDEC0,MDirection::Ref(MDirection::HADEC));
    MDirection AzEl0 = toAzEl(PhaseCenter);
    
    MDirection tmpHADec = toHADec(AzEl0);
    
    Double Az0_Rad = AzEl0.getAngle(Radian).getValue()(0);
    Double El0_Rad = AzEl0.getAngle(Radian).getValue()(1);
    
    //
    // Convert the antenna pointing offsets from (Az,El)-->(RA,Dec)-->(l,m) 
    //
    
    for(IPosition n(3,0,0,0);n(2)<=NAnt;n(2)++)
      {
	//
	// From (Az,El) -> (HA,Dec)
	//
	// Add (Az,El) offsets to the reference (Az,El)
	//
	dAz.setValue(l_off(n)+Az0_Rad);  dEl.setValue(m_off(n)+El0_Rad);
	//      dAz.setValue(0.0+Az0_Rad);  dEl.setValue(0.0+El0_Rad);
	MDirection AzEl(dAz,dEl,MDirection::Ref(MDirection::AZEL));
	//
	// Convert offsetted (Az,El) to (HA,Dec) and then to (RA,Dec)
	//
	MDirection HADec = toHADec(AzEl);
	Double HA,Dec,RA, dRA;
	HA  = HADec.getAngle(Radian).getValue()(0);
	Dec = HADec.getAngle(Radian).getValue()(1);
	RA  = LST - HA;
	dRA = RA - RA0;
	//
	// Convert offsetted (RA,Dec) -> (l,m)
	//
	l_off(n)  = sin(dRA)*cos(Dec);
	m_off(n) = sin(Dec)*CDec0-cos(Dec)*SDec0*cos(dRA);
      }
    
    return NAnt+1;
  }
  //
  //---------------------------------------------------------------
  //
  int AWProjectFT::findPointingOffsets(const VisBuffer2& vb, 
					 Cube<Float>& pointingOffsets,
					Array<Float> &l_off,
					Array<Float> &m_off,
					Bool Evaluate)
  {
    //    LogIO log_l(LogOrigin("AWProjectFT2", "findPointingOffsets[R&D]"));
    Int NAnt = 0;
    //Float tmp;
    // TBD: adapt the following to VisCal mechanism:
    MEpoch LAST;
    
    NAnt=pointingOffsets.shape()(2);
    l_off.resize(IPosition(3,2,1,NAnt));
    m_off.resize(IPosition(3,2,1,NAnt));
    IPosition ndx(3,0,0,0),ndx1(3,0,0,0);
    for(ndx(2)=0;ndx(2)<NAnt;ndx(2)++)
      {
	ndx1=ndx;
	ndx(0)=0;ndx1(0)=0;	//tmp=l_off(ndx)  = pointingOffsets(ndx1);//Axis_0,Pol_0,Ant_i
	ndx(0)=1;ndx1(0)=1;	//tmp=l_off(ndx)  = pointingOffsets(ndx1);//Axis_0,Pol_1,Ant_i
	ndx(0)=0;ndx1(0)=2;	//tmp=m_off(ndx)  = pointingOffsets(ndx1);//Axis_1,Pol_0,Ant_i
	ndx(0)=1;ndx1(0)=3;	//tmp=m_off(ndx)  = pointingOffsets(ndx1);//Axis_1,Pol_1,Ant_i
      }

    return NAnt;
    if (!Evaluate) return NAnt;
    
    //
    // Make a Coordinate Conversion Machine to go from (Az,El) to
    // (HA,Dec).
    //
    MDirection::Convert toAzEl = makeCoordinateMachine(vb,MDirection::HADEC,
						       MDirection::AZEL,
						       LAST);
    MDirection::Convert toHADec = makeCoordinateMachine(vb,MDirection::AZEL,
							MDirection::HADEC,
							LAST);
    //
    // ...and now hope that it all works and works correctly!!!
    //
    Quantity dAz(0,Radian),dEl(0,Radian);
    //
    // An array of shape [2,1,1]!
    //
    Array<Double> phaseDir = vb.subtableColumns().field().phaseDir().getColumn();
    Double RA0   = phaseDir(IPosition(3,0,0,0));
    Double Dec0  = phaseDir(IPosition(3,1,0,0));
    //  
    // Compute reference (HA,Dec)
    //
    Double LST   = LAST.get(Day).getValue();
    Double SDec0 = sin(Dec0), CDec0=cos(Dec0);
    LST -= floor(LST); // Extract the fractional day
    LST *= 2*C::pi;// Convert to Raidan
    
    Double HA0;
    HA0 = LST - RA0;
    Quantity QHA0(HA0,Radian), QDEC0(Dec0,Radian);
    //
    // Convert reference (HA,Dec) to reference (Az,El)
    //
    MDirection PhaseCenter(QHA0, QDEC0,MDirection::Ref(MDirection::HADEC));
    MDirection AzEl0 = toAzEl(PhaseCenter);
    
    MDirection tmpHADec = toHADec(AzEl0);
    
    Double Az0_Rad = AzEl0.getAngle(Radian).getValue()(0);
    Double El0_Rad = AzEl0.getAngle(Radian).getValue()(1);
    
    //
    // Convert the antenna pointing offsets from (Az,El)-->(RA,Dec)-->(l,m) 
    //
    
    for(IPosition n(3,0,0,0);n(2)<=NAnt;n(2)++)
      {
	//
	// From (Az,El) -> (HA,Dec)
	//
	// Add (Az,El) offsets to the reference (Az,El)
	//
	dAz.setValue(l_off(n)+Az0_Rad);  dEl.setValue(m_off(n)+El0_Rad);
	//      dAz.setValue(0.0+Az0_Rad);  dEl.setValue(0.0+El0_Rad);
	MDirection AzEl(dAz,dEl,MDirection::Ref(MDirection::AZEL));
	//
	// Convert offsetted (Az,El) to (HA,Dec) and then to (RA,Dec)
	//
	MDirection HADec = toHADec(AzEl);
	Double HA,Dec,RA, dRA;
	HA  = HADec.getAngle(Radian).getValue()(0);
	Dec = HADec.getAngle(Radian).getValue()(1);
	RA  = LST - HA;
	dRA = RA - RA0;
	//
	// Convert offsetted (RA,Dec) -> (l,m)
	//
	l_off(n)  = sin(dRA)*cos(Dec);
	m_off(n) = sin(Dec)*CDec0-cos(Dec)*SDec0*cos(dRA);
      }
    
    return NAnt+1;
  }
  //
  //---------------------------------------------------------------
  // The method below is for making the sensitivty image directly from
  // image-plane models of the PB.  This is efficiently done via
  // ConvolutionFunction classes, at least for the ray traced models.
  //
  // For Wide-band AWP, this is computed in a wide-band sense by
  // accumulating WT*CF and FT'ing the resulting complex grid.  This
  // is relatively more expensive but is a one-time cost.  This method
  // (defined in C++ as the name plus the function signature) is
  // therefore overloaded in AWProjectWBFT which, in case the
  // sensitivity pattern (.pb image) is not found, only issues a
  // message informing that the first gridding cycle will do this
  // accumulation and therefore will be slower than the subsequent
  // ones.  The steps necessary to compute the WB sensitivity pattern are:
  //
  //     1. Accumulate CFs (done via AWProjectWBFT::resampleCFToGrid())
  //     2. FT the grid (done via AWProjectWBFT::ftWeightImage()).
  //        This gives PBs in the feed pol basis.
  //     3. Compute the Stokes PBs from feed basis
  //
  // The accumulation is done via the
  // AWProjectWBFT::resampleCFToGrid() method (which, for accumulating
  // the CFs, is the equivalent of ::resampledDataToGrid()).  FT and
  // averaging across feed-pol planes is done in the following method:
  //
  // AWProjectWBFT::makeWBSensitivityImage(Lattice<T>&,
  //                                       ImageInterface<Float>&,
  //                                       const Matrix<Float>&,
  //                                       const Bool&)
  //
  // (which uses the AWProjectWBFT::ftWeightImage() method for FFT).
  // This method is expected to be a NoOp if avgPBReady_p state
  // variable is True.  Since these operations can be done only at the
  // _end_ of the (first) gridding cycle, they are naturally triggered
  // in the ImageInterface<Complex>& AWProjectWBFT::getImage() method,
  // which is also overloaded to first call
  // AWProjectWBFT::makeWBSensitivityImage() (which makes the PB and
  // points the class variables avgPB_p to it).  And then call
  // AWProjectFT::getImage() (which normalizes the sky images).
  //
  void AWProjectFT::makeSensitivityImage(const VisBuffer2& vb, 
					 const ImageInterface<Complex>& imageTemplate,
					 ImageInterface<Float>& sensitivityImage)
  {
    if (convFuncCtor_p->makeAverageResponse(vb, imageTemplate, sensitivityImage))
      cfCache_p->flush(sensitivityImage,sensitivityPatternQualifierStr_p); 
  }
  //
  //---------------------------------------------------------------
  //
  void AWProjectFT::makeCFPolMap(const VisBuffer2& vb, const Vector<Int>& locCfStokes,
				 Vector<Int>& polM)
  {
    Vector<Int> msStokes = vb.correlationTypes();
    Int nPol = msStokes.nelements();
    polM.resize(polMap.shape());
    polM = -1;

    for(Int i=0;i<nPol;i++)
      for(uInt j=0;j<locCfStokes.nelements();j++)
	if (locCfStokes(j) == msStokes(i))
	    {polM(i) = j;break;}
  }
  //
  //---------------------------------------------------------------
  //
  // Given a polMap (mapping of which Visibility polarization is
  // gridded onto which grid plane), make a map of the conjugate
  // planes of the grid E.g, for Stokes-I and -V imaging, the two
  // planes of the uv-grid are [LL,RR].  For input VisBuffer2
  // visibilites in order [RR,RL,LR,LL], polMap = [1,-1,-1,0].  The
  // conjugate map will be [0,-1,-1,1].
  //
  void AWProjectFT::makeConjPolMap(const VisBuffer2& vb, 
				     const Vector<Int> cfPolMap, 
				     Vector<Int>& conjPolMap)
  {
    if (conjPolMap.nelements() > 0) return;

    LogIO log_l(LogOrigin("AWProjectFT2", "makConjPolMap[R&D]"));

    //
    // All the Natak (Drama) below with slicers etc. is to extract the
    // Poln. info. for the first IF only (not much "information
    // hiding" for the code to slice arrays in a general fashion).
    //
    // Extract the shape of the array to be sliced.
    //
    //    Array<Int> stokesForAllIFs = vb.subtableColumns().polarization().corrType().getColumn();
    log_l << "############....temp code!!!!!!!!!! "
	  << SynthesisUtils::mapSpwIDToPolID(vb, vb.spectralWindows()[0])
	  << LogIO::DEBUG1;
    Vector<Int> polIDs=SynthesisUtils::mapSpwIDToPolID(vb,vb.spectralWindows()[0]);
    if (polIDs.nelements()==0)
      log_l << "Internal Error: Selected SPW did not map to any pol!" << LogIO::EXCEPTION;
    //
    // Use the first selected Spw to determine the pol. mapping in the
    // VB.  This implies that the code assumes that all selected SPWs
    // have the same pol. setup.
    //
    Int polID=polIDs(0);
    Array<Int> stokesForAllIFs = (vb.subtableColumns()).polarization().corrType().get(polID);
    IPosition stokesShape(stokesForAllIFs.shape());
    IPosition firstIFStart(stokesShape),firstIFLength(stokesShape);
    //
    // Set up the start and length IPositions to extract only the
    // first column of the array.  The following is required since the
    // array could have only one column as well.
    //
    firstIFStart(0)=0;firstIFLength(0)=stokesShape(0);
    for(uInt i=1;i<stokesShape.nelements();i++) {firstIFStart(i)=0;firstIFLength(i)=1;}
    //
    // Construct the slicer and produce the slice.  .nonDegenerate
    // required to ensure the result of slice is a pure vector.
    //
    Vector<Int> visStokes = stokesForAllIFs(Slicer(firstIFStart,firstIFLength)).nonDegenerate();

    conjPolMap = cfPolMap;
    
    Int i,j,N = cfPolMap.nelements();
    for(i=0;i<N;i++)
      if (cfPolMap[i] > -1)
	{
	if      (visStokes[i] == Stokes::RR) 
	  {
	    conjPolMap[i]=-1;
	    for(j=0;j<N;j++) if (visStokes[j] == Stokes::LL) break; 
	    conjPolMap[i]=cfPolMap[j];
	  }
	else if (visStokes[i] == Stokes::LL) 
	  {
	    conjPolMap[i]=-1;
	    for(j=0;j<N;j++) if (visStokes[j] == Stokes::RR) break; 
	    conjPolMap[i]=cfPolMap[j];
	  }
	else if (visStokes[i] == Stokes::LR) 
	  {
	    conjPolMap[i]=-1;
	    for(j=0;j<N;j++) if (visStokes[j] == Stokes::RL) break; 
	    conjPolMap[i]=cfPolMap[j];
	  }
	else if (visStokes[i] == Stokes::RL) 
	  {
	    conjPolMap[i]=-1;
	    for(j=0;j<N;j++) if (visStokes[j] == Stokes::LR) break; 
	    conjPolMap[i]=cfPolMap[j];
	  }
	}
  }
  //
  //---------------------------------------------------------------
  // Convert the CFs in the supplied CFStore to
  // A(nu)<Convolution>A(nu_*)
  //
  void AWProjectFT::makeWBCFWt(CFStore2& cfs, const Double imRefFreq)
  {
    LogIO log_l(LogOrigin("AWProjectFT2", "makeWBCFWt[R&D]"));
    log_l << "Converting WTCFs to wide-band versions" << LogIO::POST;

    Vector<Int> ant1List, ant2List;
    Vector<Quantity> paList;
    ant1List = cfs.getAnt1List();
    ant2List = cfs.getAnt2List();
    paList   = cfs.getPAList();

    if (paNdxProcessed_p.nelements() == 0) {paNdxProcessed_p.resize(1); paNdxProcessed_p[0]=false;}
    CountedPtr<CFBuffer> cfb_l, cfb_clone;
    Quantity dPA(360.0,"deg");
    for (uInt pa=0;pa<paList.nelements();pa++)
      for (uint a1=0;a1<ant1List.nelements(); a1++)
	for (uint a2=0;a2<ant2List.nelements(); a2++)
	  {
	    if (paNdxProcessed_p.nelements() < pa) {paNdxProcessed_p.resize(pa+1,true); paNdxProcessed_p[pa]=false;}
	    if (!paNdxProcessed_p[pa])
	      {
		paNdxProcessed_p[pa]=true;
		Vector<Double> wVals, fVals; PolMapType mVals, mNdx, conjMVals, conjMNdx;
		Double fIncr, wIncr;
		cfb_l = cfs.getCFBuffer(paList[pa], dPA, ant1List(a1), ant2List(a2));
		cfb_l->getCoordList(fVals,wVals,mNdx, mVals, conjMNdx, conjMVals, fIncr, wIncr);
		
		cfb_clone=cfb_l->clone();
		cfb_clone->getCoordList(fVals,wVals,mNdx, mVals, conjMNdx, conjMVals, fIncr, wIncr);

		CountedPtr<Array<Complex> > cfc_l, conjCFC_l, result_l;
		//
		// Damn!  Convolver does not work for complex convolution!
		//		Convolver<Complex> convolver;

		for (Int nw=0;nw<(Int)wVals.nelements(); nw++)
		  for (Int nf=0;nf<(Int)fVals.nelements(); nf++)
		    for (Int ipol=0;ipol<(Int)conjMVals.nelements();ipol++)
		      for (Int mRow=0;mRow<(Int)conjMVals[ipol].nelements(); mRow++)
			{
			  Double f, cf;
			  f=fVals[nf];
			  cf=sqrt(2*imRefFreq*imRefFreq - f*f);
			  Int conjNF = cfb_l->nearestFreqNdx(cf);

			  result_l  = cfb_l->getCFCellPtr(nf,nw,conjMNdx[ipol][mRow])->getStorage();
			  cfc_l     = cfb_clone->getCFCellPtr(nf,nw,conjMNdx[ipol][mRow])->getStorage();
			  conjCFC_l = cfb_clone->getCFCellPtr(conjNF,nw,conjMNdx[ipol][mRow])->getStorage();

			  SynthesisUtils::libreConvolver(*result_l,*conjCFC_l);
			}
	      }
	  }
  }
  //
  // Locate a convlution function.  It will be either in the cache
  // (mem. or disk cache) or will be computed and cached for possible
  // later use.
  //
  void AWProjectFT::findConvFunction(const ImageInterface<Complex>& image,
				     const VisBuffer2& vb)
  {
    if (!paChangeDetector.changed(vb,0)) return;
    Int cfSource=CFDefs::NOTCACHED;
    CoordinateSystem ftcoords;
    // Think of a generic call to get the key-values.  And a
    // overloadable method (or an externally supplied one?) to convert
    // the values to key-ids.  That will ensure that AWProjectFT
    // remains the A-Projection algorithm implementation configurable
    // by the behaviour of the supplied objects.
    Float pa=getVBPA(vb);
    //UUU// ok();
    visResampler_p->setMaps(chanMap, polMap); //UUU Added here.
    visResampler_p->setFreqMaps(expandedSpwFreqSel_p,expandedSpwConjFreqSel_p);

    lastPAUsedForWtImg = currentCFPA = pa;

    //Vector<Vector<Double> > pointingOffset(convFuncCtor_p->findPointingOffset(image,vb, doPointing));

    // PO::setOverSampling needs to be called here since
    // convFuncCtor_p gets that value from a combination of (1)
    // ATerm_OVERSAMPLING env. variable, (2) ATERM.OVERSAMPLING in
    // ~/.casa and (3) from existing CFCache.  This setting in the AWP
    // constructor will only get the default value from ATerm.h
    convSampling=convFuncCtor_p->getOversampling();
    po_p->setOverSampling(convFuncCtor_p->getOversampling());
    // PO::fetchPointingOffset() only updates the internal cache in PO
    // class.  PO::pullPointingOffset() is required to extract in the
    // calling class.
    po_p->fetchPointingOffset(image,vb, doPointing);

    Float dPA = paChangeDetector.getParAngleTolerance().getValue("rad");
    Quantity dPAQuant = Quantity(paChangeDetector.getParAngleTolerance());

    vb2CFBMap_p->setDoPointing(doPointing);
    cfSource = vb2CFBMap_p->makeVBRow2CFBMap(*cfs2_p,
						vb,
						dPAQuant,
						chanMap,polMap,po_p);

    if (cfSource == CFDefs::NOTCACHED)
      {
	PolMapType polMat, polIndexMat, conjPolMat, conjPolIndexMat;
	Vector<Int> visPolMap(vb.correlationTypes());
	polMat = pop_p->makePolMat(visPolMap,polMap);
	polIndexMat = pop_p->makePol2CFMat(visPolMap,polMap);

	conjPolMat = pop_p->makeConjPolMat(visPolMap,polMap);
	conjPolIndexMat = pop_p->makeConjPol2CFMat(visPolMap,polMap);

	convFuncCtor_p->setPolMap(polMap);
	convFuncCtor_p->setSpwSelection(spwChanSelFlag_p);
	convFuncCtor_p->setSpwFreqSelection(spwFreqSel_p);

	// USEFUL DEBUG MESSAGE
	//cerr << "Freq. selection: " << expandedSpwFreqSel_p << endl << expandedSpwConjFreqSel_p << endl;
	Bool pleaseDoAlsoFillTheCF=!dryRun();
	convFuncCtor_p->makeConvFunction(image,vb,wConvSize, 
					 pop_p, pa, dPA, uvScale, uvOffset,spwFreqSel_p,
					 *cfs2_p, *cfwts2_p, pleaseDoAlsoFillTheCF);
      }
    //
    // If one-time-operations in the CFCache not yet done, set the
    // pol. index maps in the CFCache.
    //
    if (!cfCache_p->OTODone())
      {
	//Vector<Int> visPolMap(vb.corrType());
	Vector<Int> visPolMap(vb.correlationTypes());

	PolMapType polMat, conjPolMat;
	polMat = pop_p->makePolMat(visPolMap,polMap);
	conjPolMat = pop_p->makeConjPolMat(visPolMap,polMap);

	PolMapType pNdx, cpNdx;
	pNdx = pop_p->makePol2CFMat(visPolMap,polMap);
	cpNdx = pop_p->makeConjPol2CFMat(visPolMap,polMap);
    
	cfCache_p->initPolMaps(pNdx,cpNdx);

	cfs2_p->initMaps(vb,spwFreqSel_p,imRefFreq_p);
	cfwts2_p->initMaps(vb,spwFreqSel_p,imRefFreq_p);
      }
    //
    // Load the average PB (sensitivity pattern) from the cache.  If
    // not found in the cache, make one and cache it.
    //
	std::tuple<int, double>cubeinfo(1,-1.0);
        double freqofBegChan;
        spectralCoord_p.toWorld(freqofBegChan, 0.0);
        
        cubeinfo=std::make_tuple(image.shape()(3),freqofBegChan);


        if(!avgPBReady_p)
          avgPBReady_p = (cfCache_p->loadAvgPB(avgPB_p,sensitivityPatternQualifierStr_p, cubeinfo) != CFDefs::NOTCACHED);
    
    if(avgPBReady_p){
        LatticeExprNode le( max( *avgPB_p ) );
        Float avgPB_max=le.getFloat();
        
        if(avgPB_max <= 0.0) avgPBReady_p = false;
    }
    
    if(!avgPBReady_p) makeSensitivityImage(vb,image,*avgPB_p);

	
    
    verifyShapes(avgPB_p->shape(), image.shape());

    if (paChangeDetector.changed(vb,0)) paChangeDetector.update(vb,0);
    //
    // Write some useful info. to the logger.
    //
    if (cfSource != CFDefs::MEMCACHE)
      {
	// If dry run, write the uvgrid as an image for later use in
	// filling the empty CFCache.  Only the co-ordinate system of
	// the uvgrid is required later.
	if (dryRun())
	  {
	    PagedImage<Complex> thisGrid(image.shape(),image.coordinates(), 
					 cfCache_p->getCacheDir()+"/uvgrid.im");
	  }
	// Save only the CF Cube for the current value of PA (not the
	// entire CFStore -- CFs for PA values encountered earlier
	// than current value have already need made persistent).
	cfs2_p->makePersistent(cfCache_p->getCacheDir().c_str(),"","",    Quantity(pa,"rad"),dPAQuant,0,0);
	cfwts2_p->makePersistent(cfCache_p->getCacheDir().c_str(),"","WT",Quantity(pa,"rad"),dPAQuant,0,0);
	Double memUsed=cfs2_p->memUsage();
	String unit(" KB");
	memUsed = (Int)(memUsed/1024.0+0.5);
	if (memUsed > 1024) {memUsed /=1024; unit=" MB";}

	LogIO log_l(LogOrigin("AWProjectFT2", "findConvFunction[R&D]"));
	log_l << "Convolution function memory footprint:" 
	      << (Int)(memUsed) << unit << " out of a maximum of "
	      << HostInfo::memoryTotal(true)/1024 << " MB" << LogIO::POST;
	
	//
	// Initialize any internal maps that may be used later for
	// efficient access.
	//
	cfs2_p->initMaps(vb,spwFreqSel_p,imRefFreq_p);
	cfwts2_p->initMaps(vb,spwFreqSel_p,imRefFreq_p);
      }
  }
  //
  //------------------------------------------------------------------------------
  // Vectorized initializeToVis.  See design related comments in the
  // .h file.
  //
  //
  // The functional goals here are:
  //
  //  1. If doPBCorrection==true (i.e. the input sky-image is flat-sky
  //  image), divide the image by the PB
  //
  //  2. Convert from Stokes to Feed (correlation) frame
  //
  //  3. FFT the image to produce gridded vis.
  //
  //  4. And since the same image buffer is used to accumulate the
  //  model, multiply the sky-image back with PB
  //
  // The call to non-vectorized version of initializeToVis() which is
  // pure virtual and hence the local version is called, is only to do
  // the FFT.  FFT is *NOT* in place.
  //
  // These operations effecitvely maintains the model image as PB*Sky
  // and supplies flat-sky model for prediction.  This can probably be
  // achieve with fewer operations and same memory buffers....but that
  // for later (SB)
  void AWProjectFT::initializeToVis(Block<CountedPtr<ImageInterface<Complex> > > & compImageVec,
				    PtrBlock<SubImage<Float> *> & modelImageVec, 
				    PtrBlock<SubImage<Float> *>& weightImageVec, 
				    PtrBlock<SubImage<Float> *>& /*fluxScaleVec*/, 
				    Block<Matrix<Float> >& weightsVec,
				    const VisBuffer2& vb)
  {
    LogIO log_p(LogOrigin("AWProjectFT2","initToVis[V][R&D]"));
    //
    // Setting the image below is crucial since init() and
    // initMaps(vb) below expect this to be set.
    //
    image=&(*compImageVec[0]);

    log_p << "Total flux in model image (before avgPB normalization): " 
	  << sum((*(modelImageVec[0])).get()) 
	  << LogIO::POST;
    if(doPBCorrection) 
      {
	// Make the sensitivity Image if applicable
	init(vb);
	initMaps(vb);
	findConvFunction(*(compImageVec[0]), vb); // Pure virtual -- call local version

	if (isDryRun) return;

	// Get the sensitivity Image
	Matrix<Float> tempWts;
	tempWts.resize();
	getWeightImage(*(weightImageVec[0]), tempWts);  // Pure virtual -- call local version

	// Normalize the model image by the sensitivity image only.
	// No local implementation -- call FTMachine version

	// Divide by avgPB ///// PBWeight
	//

	// Divide by sqrt(avgPB) ////// PBSQWeight
	//
	 normalizeImage( *(modelImageVec[0]) , weightsVec[0],  *(weightImageVec[0]), 
	 		false, (Float)pbLimit_p, (Int)4);
      }
    log_p << "Total flux in model image (after avgPB normalization): " 
	  << sum((*(modelImageVec[0])).get()) << LogIO::POST;
    
    // Convert from Stokes planes to Correlation planes
    // No local implementation -- call FTMachine version
    stokesToCorrelation(*(modelImageVec[0]), *(compImageVec[0]));

    // Call initializeToVis
    initializeToVis(*(compImageVec[0]), vb); // Pure virtual
    
    // Multiply the flat-sky model by the PB.
    // No local implementation -- call FTMachine version
    if(doPBCorrection)
      // Multiply by avgPB ///////  PBWeight
      //normalizeImage( *(modelImageVec[0]) , weightsVec[0], *(weightImageVec[0]) , false, (Float)pbLimit_p, (Int)3);
      
      // Multiply by sqrt(avgPB) ///// PBSQWeight
      normalizeImage( *(modelImageVec[0]) , weightsVec[0], *(weightImageVec[0]) , false, (Float)pbLimit_p, (Int)5);
  }
  //
  //------------------------------------------------------------------------------
  //
  void AWProjectFT::initializeToVis(ImageInterface<Complex>& iimage,
				    const VisBuffer2& vb)
  {
    LogIO log_l(LogOrigin("AWProjectFT2", "initializeToVis[R&D]"));
    image=&iimage;
    
    ok();
    
    init(vb);
    makingPSF = false;
    initMaps(vb);
    
    findConvFunction(*image, vb);
    if (isDryRun) return;
    //  
    // Initialize the maps for polarization and channel. These maps
    // translate visibility indices into image indices
    //

    nx    = image->shape()(0);
    ny    = image->shape()(1);
    npol  = image->shape()(2);
    nchan = image->shape()(3);
    
    //
    // If we are memory-based then read the image in and create an
    // ArrayLattice otherwise just use the PagedImage
    //

    isTiled=false;

      {
	IPosition gridShape(4, nx, ny, npol, nchan);
	griddedData.resize(gridShape);
	griddedData=Complex(0.0);
	
	IPosition stride(4, 1);
	IPosition blc(4, (nx-image->shape()(0)+(nx%2==0))/2,
		      (ny-image->shape()(1)+(ny%2==0))/2, 0, 0);
	IPosition trc(blc+image->shape()-stride);
	
	IPosition start(4, 0);
	griddedData(blc, trc) = image->getSlice(start, image->shape());
	
	lattice=new ArrayLattice<Complex>(griddedData);
      }

    //AlwaysAssert(lattice, AipsError);
    
    log_l << LogIO::DEBUGGING << "Starting FFT of image" << LogIO::POST;
    
    Vector<Float> sincConv(nx);
    Float centerX=nx/2;
    for (Int ix=0;ix<nx;ix++) 
      {
	Float x=C::pi*Float(ix-centerX)/(Float(nx)*Float(convSampling));
	if(ix==centerX) sincConv(ix)=1.0;
	else            sincConv(ix)=sin(x)/x;
	sincConv(ix) = 1.0;
      }
    
    if(cfCache_p->avgPBReady()) //SB
    {
      //      normalizeAvgPB();
      
      IPosition cursorShape(4, nx, 1, 1, 1);
      IPosition axisPath(4, 0, 1, 2, 3);
      LatticeStepper lsx(lattice->shape(), cursorShape, axisPath);
      LatticeIterator<Complex> lix(*lattice, lsx);
	  
      verifyShapes(avgPB_p->shape(), image->shape());
      Array<Float> avgBuf; avgPB_p->get(avgBuf);
      // If the total-power sensitivity pattern peak is too low, warn
      // the user.  This usually is indicative of a rat somewhere in
      // the pipes upstream...
      if ((sensitivityPatternQualifier_p==0) && (max(avgBuf) < 1e-04))
	log_l << "Normalization by PB requested but either PB was not"
	      <<" found in the cache or is ill-formed. Peak = "
	      << max(avgBuf)// << " " << sensitivityPatternQualifier_p
	      << LogIO::WARN << LogIO::POST;
	  

      LatticeStepper lpb(avgPB_p->shape(),cursorShape,axisPath);
      LatticeIterator<Float> lipb(*avgPB_p, lpb);

      Vector<Complex> griddedVis;
      //
      // Grid correct in anticipation of the convolution by the
      // convFunc.  Each polarization plane is corrected by the
      // appropraite primary beam.
      //
      for(lix.reset(),lipb.reset();!lix.atEnd();lix++,lipb++) 
	{
	  Int iy=lix.position()(1);
	  griddedVis = lix.rwVectorCursor();
	  
	  Vector<Float> PBCorrection(lipb.rwVectorCursor().shape());
	  PBCorrection = lipb.rwVectorCursor();
	  for(int ix=0;ix<nx;ix++) 
	    {
	      if (doPBCorrection)
		{
		  PBCorrection(ix) = pbFunc(PBCorrection(ix),pbLimit_p)*(sincConv(ix)*sincConv(iy));
		  lix.rwVectorCursor()(ix) /= (PBCorrection(ix));
		}
	      else 
		lix.rwVectorCursor()(ix) /= (1.0/(sincConv(ix)*sincConv(iy)));
	    }
	}
    }
    //
    // Now do the FFT2D in place
    //

    LatticeFFT::cfft2d(*lattice);
    log_l << LogIO::DEBUGGING << "Finished FFT" << LogIO::POST;
  }
  //
  //------------------------------------------------------------------------------
  //
  void AWProjectFT::initializeToVis(ImageInterface<Complex>& iimage,
				     const VisBuffer2& vb,
				     Array<Complex>& griddedVis,
				     Vector<Double>& uvscale)
  {
    initializeToVis(iimage, vb);
    griddedVis.assign(griddedData); //using the copy for storage
    uvscale.assign(uvScale);
    
  }
  //
  //---------------------------------------------------------------
  //
  void AWProjectFT::finalizeToVis()
  {
    visResampler_p->runTimeDG_p=0.0;
    logIO() << LogOrigin("AWProjectFT", "finalizeToVis")  << LogIO::NORMAL;
    logIO()<< LogIO::WARN << "Time degrid " << timedegrid_p << LogIO::POST;
    timedegrid_p=0.0;

  if(!lattice.null()) lattice=0;
  griddedData.resize();
  
    if(isTiled) 
      {
	AlwaysAssert(imageCache, AipsError);
	AlwaysAssert(image, AipsError);
	ostringstream o;
	imageCache->flush();
	imageCache->showCacheStatistics(o);

	LogIO log_l(LogOrigin("AWProjectFT2", "finalizeToVis[R&D]"));
	log_l << o.str() << LogIO::POST;
      }
    if(pointingToImage) delete pointingToImage;
    pointingToImage=0;
  }
  //
  //---------------------------------------------------------------
  //
  // Initialize the FFT to the Sky. Here we have to setup and
  // initialize the grid.
  //
  void AWProjectFT::initializeToSky(ImageInterface<Complex>& iimage,
				     Matrix<Float>& weight,
				     const VisBuffer2& vb)
  {
    LogIO log_l(LogOrigin("AWProjectFT2", "initializeToSky[R&D]"));
    
    // image always points to the image
    image=&iimage;
    
    init(vb);
    initMaps(vb);
    log_l << "Computed maps using FTMachine::initMaps. " << "polMap = " << polMap << LogIO::POST;
    visResampler_p->setMaps(chanMap, polMap);
    visResampler_p->setFreqMaps(expandedSpwFreqSel_p,expandedSpwConjFreqSel_p);
    
    // Initialize the maps for polarization and channel. These maps
    // translate visibility indices into image indices
    nx    = image->shape()(0);
    ny    = image->shape()(1);
    npol  = image->shape()(2);
    nchan = image->shape()(3);
    
    sumWeight=0.0;
    sumCFWeight = 0.0;
    weight.resize(sumWeight.shape());
    weight=0.0;

    //
    // Construct the HPG.  This is now donw in AWVRHPG.  Not the right
    // place, I think, but good enough for testing. (21Dec2020).
    //
    // const std::array<unsigned, 4> gridSize{(uInt)nx,(uInt)ny,(uInt)npol,(uInt)nchan};
    // const std::array<float, 2> gridScale{(float)uvScale(0), (float)uvScale(1)};
    // unsigned max_added_tasks=1; // This allocates 2 CUDA streams in AWVRHPG.
    // visResampler_p->initGridder(max_added_tasks,gridSize,gridScale);


    //
    // Initialize for in memory or to disk gridding. lattice will
    // point to the appropriate Lattice, either the ArrayLattice for
    // in memory gridding or to the image for to disk gridding.
    //
    
  
	IPosition gridShape(4, nx, ny, npol, nchan);
	if(!useDoubleGrid_p){
	    griddedData.resize(gridShape);
	  griddedData=Complex(0.0); 
	}
	else	  
  {
	  griddedData2.resize(gridShape);
	  griddedData2=DComplex(0.0);
	}
      

    //cerr << "initializeToSky for grid" << endl;
    if(useDoubleGrid_p) 
      visResampler_p->initializeToSky(griddedData2, sumWeight);
    else
      visResampler_p->initializeToSky(griddedData, sumWeight);
  }
  
  //
  //---------------------------------------------------------------
  //
  void AWProjectFT::finalizeToSky()
  {
    logIO() << LogOrigin("AWProjectFT", "finalizeToSky")  << LogIO::NORMAL;
    logIO() <<   LogIO::WARN << "time to massage data " << timemass_p << LogIO::POST;
    logIO() <<  LogIO::WARN << "time gridding " << timegrid_p << LogIO::POST;
   timemass_p=0.0;
   timegrid_p=0.0;
   if(name()=="AWProjectWBFTHPG"){
    Matrix<Double> tmpSumWgt(sumWeight.shape());
    tmpSumWgt=0.0;
    if(useDoubleGrid_p) 
      visResampler_p->finalizeToSky(griddedData2, tmpSumWgt);
    else
      visResampler_p->finalizeToSky(griddedData, tmpSumWgt);
    sumWeight=tmpSumWgt;
   
    return;
    
   }
    //
    // Now we flush the cache and report statistics For memory based,
    // we don't write anything out yet.
    //
    //    LogIO log_l(LogOrigin("AWProjectFT2", "finalizeToSky[R&D]"));

    if(pointingToImage) delete pointingToImage;
    pointingToImage=0;

    paChangeDetector.reset();
    cfCache_p->flush();
    if(useDoubleGrid_p) 
      visResampler_p->finalizeToSky(griddedData2, sumWeight);
    else
      visResampler_p->finalizeToSky(griddedData, sumWeight);
  }
  //
  //---------------------------------------------------------------
  //
  Array<Complex>* AWProjectFT::getDataPointer(const IPosition& centerLoc2D,
					       Bool readonly) 
  {
    Array<Complex>* result;
    // Is tiled: get tiles and set up offsets
    centerLoc(0)=centerLoc2D(0);
    centerLoc(1)=centerLoc2D(1);
    result=&imageCache->tile(offsetLoc, centerLoc, readonly);
    gridder->setOffset(IPosition(2, offsetLoc(0), offsetLoc(1)));
    return result;
  }
  
  // The following file has the runFORTRAN* stuff. Moving it to a
  // separate file to reduce clutter and ultimately delete it.
#include "AWProjectFT.FORTRANSTUFF"
  //
  //---------------------------------------------------------------
  //
    void AWProjectFT::put(const VisBuffer2& vb, Int /*row*/, Bool dopsf,
			FTMachine::Type type)
  {

    
    matchChannel(vb);
 

    //cerr << "CHANMAP " << chanMap << endl;
    //No point in reading data if its not matching in frequency
    if(max(chanMap)==-1)
      return;
    // Take care of translation of Bools to Integer
    makingPSF=dopsf;
    if(dopsf)
      ftmType_p=refim::FTMachine::PSF;
    Timer tim;
    tim.mark();

    
    try
      {
	findConvFunction(*image, vb);
      }
    catch(AipsError& x)
      {
	LogIO log_l(LogOrigin("AWProjectFT2", "put[R&D]"));
	log_l << x.getMesg() << LogIO::WARN;
	return;
      }
    if (isDryRun) return;

    Nant_p     = vb.subtableColumns().antenna().nrow();

    Matrix<Float> imagingweight;
    getImagingWeight(imagingweight, vb);

    Cube<Complex> data;
    //Fortran gridder need the flag as ints 
    Cube<Int> flags;
    Matrix<Float> elWeight;

    interpolateFrequencyTogrid(vb, imagingweight,data, flags, elWeight, type);

    //    
    // Get the uvws in a form that Fortran can use and do that
    // necessary phase rotation. On a Pentium Pro 200 MHz when null,
    // this step takes about 50us per uvw point. This is just barely
    // noticeable for Stokes I continuum and irrelevant for other
    // cases.
    //
    Matrix<Double> uvw(negateUV(vb));
    
    Vector<Double> dphase(vb.nRows());
    dphase=0.0;
    doUVWRotation_p=true;
    girarUVW(uvw, dphase, vb);
    refocus(uvw, vb.antenna1(), vb.antenna2(), dphase, vb);
    
    //Here we redo the match or use previous match
    
    //Channel matching for the actual spectral window of buffer

    matchChannel(vb);
    VBStore vbs;
    Vector<Int> gridShape = griddedData2.shape().asVector();
    setupVBStore(vbs,vb, elWeight,data,uvw,flags, dphase,dopsf,gridShape);
    timemass_p +=tim.real();
    tim.mark();
    if (useDoubleGrid_p)
      {
	resampleDataToGrid(griddedData2, vbs, vb, dopsf);//, *imagingweight, *data, uvw,flags,dphase,dopsf);
      }
    else
      {
	resampleDataToGrid(griddedData, vbs, vb, dopsf);//, *imagingweight, *data, uvw,flags,dphase,dopsf);
      }
    timegrid_p+=tim.real();
  }

  std::shared_ptr<std::complex<double>> AWProjectFT::getGridPtr(size_t& size) const
  {
    return visResampler_p->getGridPtr(size);
  }

  std::shared_ptr<double> AWProjectFT::getSumWeightsPtr(size_t& size) const
  {
    return visResampler_p->getSumWeightsPtr(size);
  }

  //
  //-------------------------------------------------------------------------
  // Gridding
  void AWProjectFT::resampleDataToGrid(Array<Complex>& griddedData_l, VBStore& vbs, 
				       const VisBuffer2& /*vb*/, Bool& dopsf)
  {
    visResampler_p->DataToGrid(griddedData_l, vbs, sumWeight, dopsf); 
  }
  //
  //-------------------------------------------------------------------------
  // Gridding
  void AWProjectFT::resampleDataToGrid(Array<DComplex>& griddedData_l, VBStore& vbs, 
				       const VisBuffer2& /*vb*/, Bool& dopsf)
  {
    visResampler_p->DataToGrid(griddedData_l, vbs, sumWeight, dopsf); 
  }
  //
  //---------------------------------------------------------------
  //
  void AWProjectFT::get(VisBuffer2& vb, Int /*row*/)
  {

    matchChannel(vb);
 

    //cerr << "CHANMAP " << chanMap << endl;
    //No point in reading data if its not matching in frequency
    if(max(chanMap)==-1)
      return;

    
    findConvFunction(*image, vb);
    Timer tim;
    tim.mark();
    Nant_p     = vb.subtableColumns().antenna().nrow();
    // Get the uvws in a form that Fortran can use
    Matrix<Double> uvw(negateUV(vb));

    Vector<Double> dphase(vb.nRows());
    dphase=0.0;
    doUVWRotation_p=true;
    //rotateUVW(uvw, dphase, vb);

    girarUVW(uvw, dphase, vb);
    refocus(uvw, vb.antenna1(), vb.antenna2(), dphase, vb);
    
    matchChannel(vb);
    //No point in reading data if its not matching in frequency
    if(max(chanMap)==-1) return;
    
    Cube<Complex> data;
    Cube<Int> flags;
    getInterpolateArrays(vb, data, flags);

    VBStore vbs;
    Bool tmpDoPSF=false;

    setupVBStore(vbs,vb, vb.imagingWeight(),data,uvw,flags, dphase,tmpDoPSF,griddedData.shape().asVector());

     tim.mark();
     resampleGridToData(vbs, griddedData, vb);//, uvw, flags, dphase);
     timedegrid_p+=tim.real();
    interpolateFrequencyFromgrid(vb, data, FTMachine::MODEL);
  }
  //
  //-------------------------------------------------------------------------
  // De-gridding
  void AWProjectFT::resampleGridToData(VBStore& vbs, Array<Complex>& griddedData_l,
				       const VisBuffer2& /*vb*/)
  {
    visResampler_p->GridToData(vbs, griddedData_l);
  }
  //
  //---------------------------------------------------------------
  //
  // Finalize the FFT to the Sky. Here we actually do the FFT and
  // return the resulting image
  ImageInterface<Complex>& AWProjectFT::getImage(Matrix<Float>& weights,
						  Bool fftNormalization) 
  {
    LogIO log_l(LogOrigin("AWProjectFT2", "getImage[R&D]"));

    AlwaysAssert(image, AipsError);

    weights.resize(sumWeight.shape());
    convertArray(weights, sumWeight);
    //  
    // If the weights are all zero then we cannot normalize otherwise
    // we don't care.
    //
    if(max(weights)==0.0) 
      log_l // UUU //<< LogIO::SEVERE
	    << "No useful data in " << name() << ".  Weights all zero"
	    << LogIO::POST;
    else
      {
	log_l << "Sum of weights: " << weights << " " << max(griddedData2) << " " << min(griddedData2) << LogIO::POST;
	//cerr << "Sum of weights: " << setprecision(20) << weights << endl;
      }
    // UUU else
      {
	log_l << LogIO::DEBUGGING
		<< "Starting FFT and scaling of image" << LogIO::POST;
	//    
	// x and y transforms (lattice has the gridded vis.  Make the
	// dirty images)
	//
	if (useDoubleGrid_p)
	  {
	    ArrayLattice<DComplex> darrayLattice(griddedData2);
	    LatticeFFT::cfft2d(darrayLattice,false);
	    griddedData.resize(griddedData2.shape());
	    convertArray(griddedData, griddedData2);
	    SynthesisUtilMethods::getResource("mem peak in getImage");
	
	    //Don't need the double-prec grid anymore...
	    griddedData2.resize();
	    lattice=new ArrayLattice<Complex>(griddedData);
	  }
	else
	  {
	    lattice=new ArrayLattice<Complex>(griddedData);
	    LatticeFFT::cfft2d(*lattice,false);
	  }
	const IPosition latticeShape = lattice->shape();

        int samp=getAWConvFunc()->getOversampling();
        //cerr << "SAMP " << samp << " ConvSampling "<< convSampling << endl;
        //Do sampling size correction    
        Vector<Float> sincConvX(nx);
        for (Int ix=0;ix<nx;ix++) {
          Float x=C::pi*Float(ix-nx/2)/(Float(nx)*Float(convSampling));
          if(ix==nx/2) {
            sincConvX(ix)=1.0;
          }
          else {
            sincConvX(ix)=sin(x)/x;
          }
        }
        Vector<Float> sincConvY(ny);
        for (Int ix=0;ix<ny;ix++) {
          Float x=C::pi*Float(ix-ny/2)/(Float(ny)*Float(convSampling));
          if(ix==ny/2) {
            sincConvY(ix)=1.0;
          }
          else {
            sincConvY(ix)=sin(x)/x;
          }
        }
    

        //cerr << convSampling << " max min of sincs " << max(sincConvX) << "    " << min(sincConvX) << max(sincConvY) << "     " << min(sincConvY) << endl;
	//
	// Now normalize the dirty image.
	//
	// Since *lattice is not copied to *image till the end of this
	// method, normalizeImage also needs to work with Lattices
	// (rather than ImageInterface).
	//
	// nx ny normalization from GridFT...
	{
	  Int inx = lattice->shape()(0);
	  Int iny = lattice->shape()(1);
	  Vector<Complex> correction(inx);
	  correction=Complex(1.0, 0.0);
	  // Do the Grid-correction
	  IPosition cursorShape(4, inx, 1, 1, 1);
	  IPosition axisPath(4, 0, 1, 2, 3);
	  LatticeStepper lsx(lattice->shape(), cursorShape, axisPath);
	  LatticeIterator<Complex> lix(*lattice, lsx);
	  for(lix.reset();!lix.atEnd();lix++) 
	    {
	      Int pol=lix.position()(2);
	      Int chan=lix.position()(3);
	      {

                Int iy=lix.position()(1);
                for (Int ix=0;ix<nx;ix++) {
                  correction(ix)=1.0/(sincConvX(ix)*sincConvY(iy));
                 }
                //cerr << iy << " min max corr " << min(abs(correction)) << "    " << max(abs(correction)) << endl;
                lix.rwVectorCursor()*=correction;
		if(fftNormalization) 
		  {
		    if(weights(pol,chan)!=0.0)
		      {
			Complex rnorm(Float(inx)*Float(iny)/( weights(pol,chan) ));
			lix.rwCursor()*=rnorm;
		      }
		    else
		      lix.woCursor()=0.0;
		  }
		else 
		  {
		    Complex rnorm(Float(inx)*Float(iny));
		    lix.rwCursor()*=rnorm;
		  }
	      }
	    }
	}
	if(!isTiled) 
	  {
	    //
	    // Check the section from the image BEFORE converting to a lattice 
	    //
            LatticeLocker lock1 (*(image), FileLocker::Write);
	    IPosition blc(4, (nx-image->shape()(0)+(nx%2==0))/2,
			  (ny-image->shape()(1)+(ny%2==0))/2, 0, 0);
	    IPosition stride(4, 1);
	    IPosition trc(blc+image->shape()-stride);
	    //
	    // Do the copy
	    //
	    image->put(griddedData(blc, trc));
	    

	    if(!lattice.null()) lattice=0;
	    griddedData.resize(IPosition(1,0));
	  }
      }

    return *image;
  }
  //
  //---------------------------------------------------------------
  //
  // Get weight image
  void AWProjectFT::getWeightImage(ImageInterface<Float>& weightImage,
				   Matrix<Float>& weights) 
  {
    weights.resize(sumWeight.shape());
    convertArray(weights,sumWeight);
    
    const IPosition latticeShape = weightImage.shape();
    const IPosition avgpbShape = avgPB_p->shape();

    //cout << "AWP::getWeightImage : weightimage shape : " << latticeShape << "  and avgpb shape : " << avgpbShape << " nelems " << avgpbShape.nelements()<< "  " << sumWeight << endl;
     if(avgpbShape.nelements()==0 || ( avgpbShape != latticeShape) )
      avgPB_p->resize(weightImage.shape());
    
    Int nx=latticeShape(0);
    Int ny=latticeShape(1);

    int samp=getAWConvFunc()->getOversampling();
    //cerr << "2 samp " << samp << " convSamp " << convSampling << endl;
    //Do sampling size correction    
    Vector<Float> sincConvX(nx);
    for (Int ix=0;ix<nx;ix++) {
      Float x=C::pi*Float(ix-nx/2)/(Float(nx)*Float(convSampling));
      if(ix==nx/2) {
        sincConvX(ix)=1.0;
      }
      else {
        sincConvX(ix)=sin(x)/x;
      }
    }
    Vector<Float> sincConvY(ny);
    for (Int ix=0;ix<ny;ix++) {
      Float x=C::pi*Float(ix-ny/2)/(Float(ny)*Float(convSampling));
      if(ix==ny/2) {
        sincConvY(ix)=1.0;
      }
      else {
        sincConvY(ix)=sin(x)/x;
      }
    }
    

       

    {
      IPosition cursorShape(4, nx, ny, latticeShape(2), latticeShape(3));
      IPosition axisPath(4, 0, 1, 2, 3);
      LatticeStepper lsx(latticeShape, cursorShape, axisPath);
      LatticeIterator<Float> lix(weightImage, lsx);
      LatticeIterator<Float> liy(*avgPB_p,lsx);
      for(lix.reset();!lix.atEnd();lix++) 
        {
          lix.rwCursor()=liy.cursor();
        }
    }
    {//sampling size correction
      Vector<Float> correction(nx);
      correction=1.0;
      // Do the Grid-correction
      IPosition cursorShape(4, nx, 1, 1, 1);
      IPosition axisPath(4, 0, 1, 2, 3);
      LatticeStepper lsx(weightImage.shape(), cursorShape, axisPath);
      LatticeIterator<Float> lix(weightImage, lsx);
      for(lix.reset();!lix.atEnd();lix++) 
        {
               
          Int iy=lix.position()(1);
          for (Int ix=0;ix<nx;ix++) {
            correction(ix)=1.0/(sincConvX(ix)*sincConvY(iy));
          }
          lix.rwVectorCursor()*=correction;
        }
        }
  }
  //---------------------------------------------------------------
    void AWProjectFT::setWeightImage(ImageInterface<Float>& weightImage){
      //cerr <<"@@@loading weightimage" << endl;
      IPosition latticeShape = weightImage.shape();
      CoordinateSystem cs=weightImage.coordinates();
      avgPB_p=new TempImage<Float>(latticeShape, cs);
      avgPB_p->copyData(weightImage);
      avgPBReady_p=True;

    }
    
  //---------------------------------------------------------------
  //
  Bool AWProjectFT::toRecord(RecordInterface& outRec, Bool withImage) 
  {
    // Save the current AWProjectFT object to an output state record
    Bool retval = true;
    String error;
    //save the base class variables
    if(!FTMachine::toRecord(error, outRec, withImage, ""))
      return false;
    Double cacheVal=(Double) cachesize;
    outRec.define("cache", cacheVal);
    outRec.define("tile", tilesize);
    
    Vector<Double> phaseValue(2);
    String phaseUnit;
    phaseValue=mTangent_p.getAngle().getValue();
    phaseUnit= mTangent_p.getAngle().getUnit();
    outRec.define("phasevalue", phaseValue);
    outRec.define("phaseunit", phaseUnit);
    
    Vector<Double> dirValue(3);
    String dirUnit;
    dirValue=mLocation_p.get("m").getValue();
    dirUnit=mLocation_p.get("m").getUnit();
    outRec.define("dirvalue", dirValue);
    outRec.define("dirunit", dirUnit);
    
    outRec.define("padding", padding_p);
    outRec.define("maxdataval", maxAbsData);
    
    Vector<Int> center_loc(4), offset_loc(4);
    for (Int k=0; k<4 ; k++)
      {
	center_loc(k)=centerLoc(k);
	offset_loc(k)=offsetLoc(k);
      }
    outRec.define("centerloc", center_loc);
    outRec.define("offsetloc", offset_loc);
    outRec.define("sumofweights", sumWeight);
    outRec.define("sumofcfweights", sumCFWeight);
    if(withImage && image)
      { 
	ImageInterface<Complex>& tempimage(*image);
	Record imageContainer;
	String error;
	retval = (retval || tempimage.toRecord(error, imageContainer));
	outRec.defineRecord("image", imageContainer);
      }
    return retval;
  }
  //
  //---------------------------------------------------------------
  //
  Bool AWProjectFT::fromRecord(const RecordInterface& inRec)
  {
    Bool retval = true;
    String error;
    if(!FTMachine::fromRecord(error, inRec))
      return false;
    imageCache=0; lattice=0; 
    Double cacheVal;
    inRec.get("cache", cacheVal);
    cachesize=(Long) cacheVal;
    inRec.get("tile", tilesize);
    
    Vector<Double> phaseValue(2);
    inRec.get("phasevalue",phaseValue);
    String phaseUnit;
    inRec.get("phaseunit",phaseUnit);
    Quantity val1(phaseValue(0), phaseUnit);
    Quantity val2(phaseValue(1), phaseUnit); 
    MDirection phasecenter(val1, val2);
    
    mTangent_p=phasecenter;
    // This should be passed down too but the tangent plane is 
    // expected to be specified in all meaningful cases.
    tangentSpecified_p=true;  
    Vector<Double> dirValue(3);
    String dirUnit;
    inRec.get("dirvalue", dirValue);
    inRec.get("dirunit", dirUnit);
    MVPosition dummyMVPos(dirValue(0), dirValue(1), dirValue(2));
    MPosition mLocation(dummyMVPos, MPosition::ITRF);
    mLocation_p=mLocation;
    
    inRec.get("padding", padding_p);
    inRec.get("maxdataval", maxAbsData);
    
    Vector<Int> center_loc(4), offset_loc(4);
    inRec.get("centerloc", center_loc);
    inRec.get("offsetloc", offset_loc);
    uInt ndim4 = 4;
    centerLoc=IPosition(ndim4, center_loc(0), center_loc(1), center_loc(2), 
			center_loc(3));
    offsetLoc=IPosition(ndim4, offset_loc(0), offset_loc(1), offset_loc(2), 
			offset_loc(3));
    inRec.get("sumofweights", sumWeight);
    inRec.get("sumofcfweights", sumCFWeight);
    if(inRec.nfields() > 12 )
      {
	Record imageAsRec=inRec.asRecord("image");
	if(!image) image= new TempImage<Complex>(); 
	
	String error;
	retval = (retval || image->fromRecord(error, imageAsRec));    
	
	// Might be changing the shape of sumWeight
	//init(vb); 
	
	if(isTiled) 
    	  lattice=CountedPtr<Lattice<Complex> > (image, false);
	else 
	  {
	    //
	    // Make the grid the correct shape and turn it into an
	    // array lattice Check the section from the image BEFORE
	    // converting to a lattice
	    //
	    IPosition gridShape(4, nx, ny, npol, nchan);
	    griddedData.resize(gridShape);
	    griddedData=Complex(0.0);
	    IPosition blc(4, (nx-image->shape()(0)+(nx%2==0))/2,
			  (ny-image->shape()(1)+(ny%2==0))/2, 0, 0);
	    IPosition start(4, 0);
	    IPosition stride(4, 1);
	    IPosition trc(blc+image->shape()-stride);
	    griddedData(blc, trc) = image->getSlice(start, image->shape());
	    
	    lattice=new ArrayLattice<Complex>(griddedData);
	  }
	
	AlwaysAssert(image, AipsError);
      };
    return retval;
  }
  //
  //---------------------------------------------------------------
  //
  void AWProjectFT::ok() 
  {
    AlwaysAssert(image, AipsError);
  }
  //
  //-------------------------------------------------------------------------
  //  
  void AWProjectFT::setPAIncrement(const Quantity& computePAIncrement,
				   const Quantity& rotateOTFPAIncrement)
  {
    LogIO log_l(LogOrigin("AWProjectFT2", "setPAIncrement[R&D]"));

    rotateOTFPAIncr_p = rotateOTFPAIncrement.getValue("rad");
    computePAIncr_p = computePAIncrement.getValue("rad");
    convFuncCtor_p->setRotateCF(computePAIncr_p, rotateOTFPAIncr_p);

    paChangeDetector.setTolerance(computePAIncrement);
    visResampler_p->setPATolerance(computePAIncrement.getValue("rad"));
    log_l << LogIO::NORMAL <<"Setting PA increment to " 
	  << computePAIncrement.getValue("deg") << " deg" << endl;
    cfCache_p->setPAChangeDetector(paChangeDetector);
  }
  //
  //-------------------------------------------------------------------------
  //  
  Bool AWProjectFT::verifyShapes(IPosition pbShape, IPosition skyShape)
  {
    LogIO log_l(LogOrigin("AWProjectFT2", "verifyShapes[R&D]"));

    if ((pbShape(0) != skyShape(0)) && // X-axis
	(pbShape(1) != skyShape(1)) && // Y-axis
	(pbShape(2) != skyShape(2)))   // Poln-axis
      {
	log_l << "Sky and/or polarization shape of the avgPB and the sky model do not match."
	      << LogIO::EXCEPTION;
	return false;
      }
    return true;
    
  }
  //
  //-------------------------------------------------------------------------
  //  
  void AWProjectFT::setupVBStore(VBStore& vbs,
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
    makeCFPolMap(vb,cfStokes_p,CFMap_p);
    makeConjPolMap(vb,CFMap_p,ConjCFMap_p);

    visResampler_p->setParams(uvScale,uvOffset,dphase);
    visResampler_p->setMaps(chanMap, polMap);
    visResampler_p->setCFMaps(CFMap_p, ConjCFMap_p);
    visResampler_p->setFreqMaps(expandedSpwFreqSel_p,expandedSpwConjFreqSel_p);
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
    vbs.paQuant_p = Quantity(getPA(vb),"rad");

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

    //timer_p.mark();

    po_p->fetchPointingOffset(*image, vb, doPointing);
    
    // Run the garbage collector for the supplied CFStore2 and make VB to CFB map in VB2CFBMap.
    auto cleanup_setup = [&](CountedPtr<CFStore2>& cfs_l, const VisBuffer2& vb_l)
		     {
		       cfs_l->invokeGC(vbs.spwID_p);
		       vb2CFBMap_p->setDoPointing(doPointing);
		       vb2CFBMap_p->makeVBRow2CFBMap(*cfs_l,
						     vb_l,
						     paChangeDetector.getParAngleTolerance(),
						     chanMap,polMap,po_p);
		       
		     };
    if (makingPSF || (vbs.ftmType_p==casa::refim::FTMachine::WEIGHT))
      cleanup_setup(cfwts2_p, vb);
    else
      {
	// If the Wt. CFs are still in the memory, clear them.  They
	// won't be required again (though with the silly check below,
	// if the in-memory Wt. CFs are less than 1KB, they will be
	// left in memory).
	if (cfwts2_p->memUsage() > 1000) cfwts2_p->clear();

	cleanup_setup(cfs2_p,vb);
      }

    //
    // For AzElApertures, this rotates the CFs.
    //
    convFuncCtor_p->prepareConvFunction(vb,*vb2CFBMap_p);
    
    vbs.accumCFs_p=((vbs.uvw_p.nelements() == 0) && dopsf);
    visResampler_p->setVB2CFMap(vb2CFBMap_p);
    
    // The following code is required only for GPU or multi-threaded
    //gridder.  Currently does not work without the rest of the
    //GPU/multi-threaded infrastructure (though, I (SB) thought this
    //was designed to be benign for normal gridding).
    //
    visResampler_p->initializeDataBuffers(vbs);
  }
  //
  //---------------------------------------------------------------
  //
  void AWProjectFT::get(VisBuffer2& vb, Cube<Complex>& modelVis, 
			 Array<Complex>& griddedVis, Vector<Double>& scale,
			 Int row)
  {

    (void)scale; //Suppress the warning

    Int nX=griddedVis.shape()(0);
    Int nY=griddedVis.shape()(1);
    Vector<Double> offset(2);
    offset(0)=Double(nX)/2.0;
    offset(1)=Double(nY)/2.0;
    // If row is -1 then we pass through all rows
    Int startRow, endRow, nRow;
    if (row==-1) 
      {
	nRow=vb.nRows();
	startRow=0;
	endRow=nRow-1;
	modelVis.set(Complex(0.0,0.0));
      } 
    else 
      {
	nRow=1;
	startRow=row;
	endRow=row;
      }
    
    Int NAnt=0;
    
    if (doPointing) 
      NAnt = findPointingOffsets(vb,l_offsets, m_offsets,true);
    
    
    //  
    // Get the uvws in a form that Fortran can use
    //
    Matrix<Double> uvw(negateUV(vb));
    Vector<Double> dphase(vb.nRows());
    dphase=0.0;
    doUVWRotation_p=true;
    //rotateUVW(uvw, dphase, vb);
    girarUVW(uvw, dphase, vb); 
    refocus(uvw, vb.antenna1(), vb.antenna2(), dphase, vb);
    
    // This is the convention for dphase
    //dphase*=-1.0;
    
    Cube<Int> flags(vb.flagCube().shape());
    flags=0;
    flags(vb.flagCube())=true;
    
    
    matchChannel(vb);
    Vector<Int> rowFlags(vb.nRows());
    rowFlags=0;
    rowFlags(vb.flagRow())=true;
    if(!usezero_p) 
      for (Int rownr=startRow; rownr<=endRow; rownr++) 
	if(vb.antenna1()(rownr)==vb.antenna2()(rownr)) rowFlags(rownr)=1;
    
    visResampler_p->setParams(uvScale,uvOffset,dphase);
    visResampler_p->setMaps(chanMap, polMap);
    visResampler_p->setFreqMaps(expandedSpwFreqSel_p,expandedSpwConjFreqSel_p);

    IPosition s(modelVis.shape());
    Int Conj=0,doGrad=0,ScanNo=0;
    Double area=1.0;
    Int tmpPAI=1;
    Cube<Complex> visCubeModel(vb.visCubeModel());
    runFortranGet(uvw,dphase,visCubeModel,s,Conj,flags,rowFlags,row,
		  offset,&griddedVis,nx,ny,npol,nchan,vb,NAnt,ScanNo,sigma,
		  l_offsets,m_offsets,area,doGrad,tmpPAI);
    vb.setVisCubeModel(visCubeModel);
  }
  //
  //----------------------------------------------------------------------
  //
  void AWProjectFT::initVisBuffer(VisBuffer2& vb, Type whichVBColumn)
  {
    if (whichVBColumn      == FTMachine::MODEL)    vb.setVisCubeModel(Complex(0.0,0.0));
    else if (whichVBColumn == FTMachine::OBSERVED) vb.setVisCube(Complex(0.0,0.0));
  }
  //
  //----------------------------------------------------------------------
  //
    void AWProjectFT::initVisBuffer(VisBuffer2&, // vb
				    Type, // whichVBColumn
				    Int // row
				    )
  {
    cerr << "AWProjectFT::initVisBuffer(VisBuffer2& vb, Type whichVBColumn, Int row) disabled" << endl;
  }

    void AWProjectFT::ComputeResiduals(VisBuffer2&vb, Bool useCorrected)
    {
      VBStore vbs;

      vbs.nRow_p = vb.nRows();
      vbs.beginRow_p = 0;
      vbs.endRow_p = vbs.nRow_p;

      vbs.modelCube_p.reference(vb.visCubeModel());
      if (useCorrected) vbs.correctedCube_p.reference(vb.visCubeCorrected());
      else vbs.visCube_p.reference(vb.visCube());
      vbs.useCorrected_p = useCorrected;
      visResampler_p->ComputeResiduals(vbs);
    }

    //
    //---------------------------------------------------------------
    //---------------------------------------------------------------
    //---------------------------------------------------------------
    //---------------------------------------------------------------
    // THIS IS FOR NON-PRODUCTION WORK.  HERE SINCE IT DOES LINK TO POINTING SELFCAL CODE
    //
    // Predict the coherences as well as their derivatives w.r.t. the
    // pointing offsets.
    //
    void AWProjectFT::nget(VisBuffer2& vb,
			   // These offsets should be appropriate for the VB
			   Array<Float>& l_off, Array<Float>& m_off,
			   Cube<Complex>& Mout,
			   Cube<Complex>& dMout1,
			   Cube<Complex>& dMout2,
			   Int Conj, Int doGrad)
    {
      LogIO log_l(LogOrigin("AWProjectFT2", "nget[R&D]"));
      Int startRow, endRow, nRow;
      nRow=vb.nRows();
      startRow=0;
      endRow=nRow-1;

      Mout = dMout1 = dMout2 = Complex(0,0);

      findConvFunction(*image, vb);
      Int NAnt=0;
      Nant_p     = vb.subtableColumns().antenna().nrow();
      if (doPointing)   
	NAnt = findPointingOffsets(vb,l_offsets,m_offsets,false);

      l_offsets=l_off;
      m_offsets=m_off;
      Matrix<Double> uvw(negateUV(vb));
    
      Vector<Double> dphase(vb.nRows());
      dphase=0.0;
      doUVWRotation_p=true;
      //rotateUVW(uvw, dphase, vb);
      girarUVW(uvw, dphase, vb);
      refocus(uvw, vb.antenna1(), vb.antenna2(), dphase, vb);
      // This is the convention for dphase
      //    dphase*=-1.0;

      Cube<Int> flags(vb.flagCube().shape());
      flags=0;
      flags(vb.flagCube())=true;
    
      Vector<Int> rowFlags(vb.nRows());
      rowFlags=0;
      rowFlags(vb.flagRow())=true;
      if(!usezero_p) 
	for (Int rownr=startRow; rownr<=endRow; rownr++) 
	  if(vb.antenna1()(rownr)==vb.antenna2()(rownr)) rowFlags(rownr)=1;
    
      IPosition s,gradS;
      Cube<Complex> visdata,gradVisAzData,gradVisElData;
      //
      // visdata now references the Mout data structure rather than to the internal VB storeage.
      //
      visdata.reference(Mout);

      if (doGrad)
	{
	  // The following should reference some slice of dMout?
	  gradVisAzData.reference(dMout1);
	  gradVisElData.reference(dMout2);
	}
      visResampler_p->setParams(uvScale,uvOffset,dphase);
      visResampler_p->setMaps(chanMap, polMap);
      visResampler_p->setFreqMaps(expandedSpwFreqSel_p,expandedSpwConjFreqSel_p);

      makeCFPolMap(vb,cfStokes_p,CFMap_p);
      makeConjPolMap(vb,CFMap_p,ConjCFMap_p);
      visResampler_p->setCFMaps(CFMap_p, ConjCFMap_p);
      //
      // Begin the actual de-gridding.
      //
      if(isTiled) 
	{
	  log_l << "The sky model is tiled" << LogIO::NORMAL << LogIO::POST;
	  Double invLambdaC=vb.getFrequencies(0)(0)/C::c;
	  Vector<Double> uvLambda(2);
	  Vector<Int> centerLoc2D(2);
	  centerLoc2D=0;
	
	  // Loop over all rows
	  for (Int rownr=startRow; rownr<=endRow; rownr++) 
	    {
	    
	      // Calculate uvw for this row at the center frequency
	      uvLambda(0)=uvw(0, rownr)*invLambdaC;
	      uvLambda(1)=uvw(1, rownr)*invLambdaC;
	      centerLoc2D=gridder->location(centerLoc2D, uvLambda);
	    
	      // Is this point on the grid?
	      if(gridder->onGrid(centerLoc2D)) 
		{
		
		  // Get the tile
		  Array<Complex>* dataPtr=getDataPointer(centerLoc2D, true);
		  gridder->setOffset(IPosition(2, offsetLoc(0), offsetLoc(1)));
		  Int aNx=dataPtr->shape()(0);
		  Int aNy=dataPtr->shape()(1);
		
		  // Now use FORTRAN to do the gridding. Remember to 
		  // ensure that the shape and offsets of the tile are 
		  // accounted for.
		
		  Vector<Double> actualOffset(3);
		  for (Int i=0;i<2;i++) 
		    actualOffset(i)=uvOffset(i)-Double(offsetLoc(i));
		
		  actualOffset(2)=uvOffset(2);
		  IPosition s(vb.visCubeModel().shape());
		
		  Int ScanNo=0, tmpPAI;
		  Double area=1.0;
		  tmpPAI = 1;
		  runFortranGetGrad(uvw,dphase,visdata,s,
				    gradVisAzData,gradVisElData,
				    Conj,flags,rowFlags,rownr,
				    actualOffset,dataPtr,aNx,aNy,npol,nchan,vb,NAnt,ScanNo,sigma,
				    l_offsets,m_offsets,area,doGrad,tmpPAI);
		}
	    }
	}
      else 
	{
	  IPosition s(vb.visCubeModel().shape());
	  Int ScanNo=0, tmpPAI, trow=-1;
	  Double area=1.0;
	  tmpPAI = 1;
	  runFortranGetGrad(uvw,dphase,visdata/*vb.modelVisCube()*/,s,
			    gradVisAzData, gradVisElData,
			    Conj,flags,rowFlags,trow,
			    uvOffset,&griddedData,nx,ny,npol,nchan,vb,NAnt,ScanNo,sigma,
			    l_offsets,m_offsets,area,doGrad,tmpPAI);
	}
    
    }
    void AWProjectFT::get(VisBuffer2& vb,       
			  VisBuffer2& gradVBAz,
			  VisBuffer2& gradVBEl,
			  Cube<Float>& pointingOffsets,
			  Int row,  // default row=-1 
			  Type whichVBColumn, // default whichVBColumn = FTMachine::MODEL
			  Type whichGradVBColumn,// default whichGradVBColumn = FTMachine::MODEL
			  Int Conj, Int doGrad) // default Conj=0, doGrad=1
    {
      // If row is -1 then we pass through all rows
      Int startRow, endRow, nRow;
      if (row==-1) 
	{
	  nRow=vb.nRows();
	  startRow=0;
	  endRow=nRow-1;
	  initVisBuffer(vb,whichVBColumn);
	  if (doGrad)
	    {
	      initVisBuffer(gradVBAz, whichGradVBColumn);
	      initVisBuffer(gradVBEl, whichGradVBColumn);
	    }
	}
      else 
	{
	  nRow=1;
	  startRow=row;
	  endRow=row;
	  initVisBuffer(vb,whichVBColumn,row);
	  if (doGrad)
	    {
	      initVisBuffer(gradVBAz, whichGradVBColumn,row);
	      initVisBuffer(gradVBEl, whichGradVBColumn,row);
	    }
	}
    
      findConvFunction(*image, vb);

      Nant_p     = vb.subtableColumns().antenna().nrow();
      Int NAnt=0;
      if (doPointing)   
	NAnt = findPointingOffsets(vb,pointingOffsets,l_offsets,m_offsets,false);

      Matrix<Double> uvw(negateUV(vb));
    
      Vector<Double> dphase(vb.nRows());
      dphase=0.0;
      doUVWRotation_p=true;
      // rotateUVW(uvw, dphase, vb);
      girarUVW(uvw, dphase, vb);
      refocus(uvw, vb.antenna1(), vb.antenna2(), dphase, vb);
    
      // This is the convention for dphase
      //dphase*=-1.0;
    
    
      Cube<Int> flags(vb.flagCube().shape());
      flags=0;
      flags(vb.flagCube())=true;
      //    
      //Check if ms has changed then cache new spw and chan selection
      //
      // if(vb.newMS()) matchAllSpwChans(vb);
    
      //Here we redo the match or use previous match
      //
      //Channel matching for the actual spectral window of buffer
      //
      // if(doConversion_p[vb.spectralWindow()])
      //   matchChannel(vb.spectralWindow(), vb);
      // else
      //   {
      // 	chanMap.resize();
      // 	chanMap=multiChanMap_p[vb.spectralWindow()];
      //   }
    
      matchChannel(vb);
      Vector<Int> rowFlags(vb.nRows());
      rowFlags=0;
      rowFlags(vb.flagRow())=true;
      if(!usezero_p) 
	for (Int rownr=startRow; rownr<=endRow; rownr++) 
	  if(vb.antenna1()(rownr)==vb.antenna2()(rownr)) rowFlags(rownr)=1;
	
      for (Int rownr=startRow; rownr<=endRow; rownr++) 
	if (vb.antenna1()(rownr) != vb.antenna2()(rownr)) 
	  rowFlags(rownr) = (vb.flagRow()(rownr)==true);
    
      IPosition s,gradS;
      Cube<Complex> visdata,gradVisAzData,gradVisElData;
      if (whichVBColumn == FTMachine::MODEL) 
	{
	  s = vb.visCubeModel().shape();
	  visdata.reference(vb.visCubeModel());
	}
      else if (whichVBColumn == FTMachine::OBSERVED)  
	{
	  s = vb.visCube().shape();
	  visdata.reference(vb.visCube());
	}
    
      if (doGrad)
	{
	  if (whichGradVBColumn == FTMachine::MODEL) 
	    {
	      //	    gradS = gradVBAz.modelVisCube().shape();
	      gradVisAzData.reference(gradVBAz.visCubeModel());
	      gradVisElData.reference(gradVBEl.visCubeModel());
	    }
	  else if (whichGradVBColumn == FTMachine::OBSERVED)  
	    {
	      //	    gradS = gradVBAz.visCube().shape();
	      gradVisAzData.reference(gradVBAz.visCube());
	      gradVisElData.reference(gradVBEl.visCube());
	    }
	}
      visResampler_p->setParams(uvScale,uvOffset,dphase);
      visResampler_p->setFreqMaps(expandedSpwFreqSel_p,expandedSpwConjFreqSel_p);
      visResampler_p->setMaps(chanMap, polMap);
      //  Vector<Int> ConjCFMap, CFMap;
      makeCFPolMap(vb,cfStokes_p,CFMap_p);
      makeConjPolMap(vb,CFMap_p,ConjCFMap_p);
      visResampler_p->setCFMaps(CFMap_p, ConjCFMap_p);
    
      if(isTiled) 
	{
	  Double invLambdaC=vb.getFrequencies(0)(0)/C::c;
	  Vector<Double> uvLambda(2);
	  Vector<Int> centerLoc2D(2);
	  centerLoc2D=0;
	
	  // Loop over all rows
	  for (Int rownr=startRow; rownr<=endRow; rownr++) 
	    {
	    
	      // Calculate uvw for this row at the center frequency
	      uvLambda(0)=uvw(0, rownr)*invLambdaC;
	      uvLambda(1)=uvw(1, rownr)*invLambdaC;
	      centerLoc2D=gridder->location(centerLoc2D, uvLambda);
	    
	      // Is this point on the grid?
	      if(gridder->onGrid(centerLoc2D)) 
		{
		
		  // Get the tile
		  Array<Complex>* dataPtr=getDataPointer(centerLoc2D, true);
		  gridder->setOffset(IPosition(2, offsetLoc(0), offsetLoc(1)));
		  Int aNx=dataPtr->shape()(0);
		  Int aNy=dataPtr->shape()(1);
		
		  // Now use FORTRAN to do the gridding. Remember to 
		  // ensure that the shape and offsets of the tile are 
		  // accounted for.
		
		  Vector<Double> actualOffset(3);
		  for (Int i=0;i<2;i++) 
		    actualOffset(i)=uvOffset(i)-Double(offsetLoc(i));
		
		  actualOffset(2)=uvOffset(2);
		  IPosition s(vb.visCubeModel().shape());
		
		  Int ScanNo=0, tmpPAI;
		  Double area=1.0;
		  tmpPAI = 1;
		  runFortranGetGrad(uvw,dphase,visdata,s,
				    gradVisAzData,gradVisElData,
				    Conj,flags,rowFlags,rownr,
				    actualOffset,dataPtr,aNx,aNy,npol,nchan,vb,NAnt,ScanNo,sigma,
				    l_offsets,m_offsets,area,doGrad,tmpPAI);
		}
	    }
	}
      else 
	{
	
	  IPosition s(vb.visCubeModel().shape());
	  Int ScanNo=0, tmpPAI;
	  Double area=1.0;

	  tmpPAI = 1;

	  runFortranGetGrad(uvw,dphase,visdata/*vb.modelVisCube()*/,s,
			    gradVisAzData, gradVisElData,
			    Conj,flags,rowFlags,row,
			    uvOffset,&griddedData,nx,ny,npol,nchan,vb,NAnt,ScanNo,sigma,
			    l_offsets,m_offsets,area,doGrad,tmpPAI);
	}
    }
    //---------------------------------------------------------------
    //---------------------------------------------------------------
    //---------------------------------------------------------------
    //---------------------------------------------------------------
    // THIS IS FOR NON-PRODUCTION WORK.  HERE SINCE IT DOES LINK TO POINTING SELFCAL CODE
    
} //# NAMESPACE CASA - END
};
