//# FringeJones.h: Declaration of fringe-fitting VisCal
//# Copyright (C) 1996,1997,2000,2001,2002,2003,2011,2016
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
//# Correspondence concerning AIPS++ should be adressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//#

#ifndef SYNTHESIS_FRINGEJONES_H
#define SYNTHESIS_FRINGEJONES_H

#include <casa/aips.h>
#include <synthesis/MeasurementComponents/StandardVisCal.h>
#include <synthesis/MeasurementComponents/DelayRateFFT.h>
#include <synthesis/CalTables/CTTimeInterp1.h>


namespace casa { //# NAMESPACE CASA - BEGIN

// Rate-aware time interpolation engine
class CTRateAwareTimeInterp1 : public CTTimeInterp1 {

public:

  // From NewCalTable
  CTRateAwareTimeInterp1(NewCalTable& ct,
			 const casacore::String& timetype,
			 casacore::Array<casacore::Float>& result,
			 casacore::Array<casacore::Bool>& rflag);

  // Destructor
  virtual ~CTRateAwareTimeInterp1();

  // Interpolate, given timestamp; returns T if new result
  virtual casacore::Bool interpolate(casacore::Double time);

  // static factory method to make CTRateAwareTimeInterp1
  // NB: returns pointer to a _generic_ CTTimeInterp1
  static CTTimeInterp1* factory(NewCalTable& ct,
				const casacore::String& timetype,
				casacore::Array<casacore::Float>& result,
				casacore::Array<casacore::Bool>& rflag) {
    return new casa::CTRateAwareTimeInterp1(ct,timetype,result,rflag); }

private:

  // Refine time-dep phase with rate info
  void applyPhaseRate(casacore::Bool single);

  // Default ctor is unusable
  //  CTRateAwareTimeInterp1(); 

};


// Fringe-fitting (parametrized phase) VisCal
class FringeJones : public GJones {
public:
  // Constructor
  //  TBD:  MSMetaInfoForCal-aware version; deprecate older ones
  FringeJones(VisSet& vs);
  FringeJones(casacore::String msname,casacore::Int MSnAnt,casacore::Int MSnSpw);
  FringeJones(const MSMetaInfoForCal& msmc);
  FringeJones(casacore::Int nAnt);

  virtual ~FringeJones();

  // We have casacore::Float parameters
  virtual VisCalEnum::VCParType parType() { return VisCalEnum::REAL; };

  // Return the type enum
  virtual Type type() { return VisCal::K; };

  // Return type name as string
  virtual casacore::String typeName()     { return "Fringe Jones"; };
  virtual casacore::String longTypeName() { return "Fringe Jones (parametrized phase)"; };

  // Type of Jones matrix according to nPar()
  virtual Jones::JonesType jonesType() { return Jones::Diagonal; };

  virtual bool timeDepMat() { return true; };

  // Freq dependence (delays)
  virtual casacore::Bool freqDepPar() { return false; };
  // If the following is false frequency spectrum is squashed to one point.
  virtual casacore::Bool freqDepMat() { return true; };

  // Local setApply to enforce calWt=F for delays
  virtual void setApply(const casacore::Record& apply);
  virtual void setApply();
  using GJones::setApply;
  virtual void setCallib(const casacore::Record& callib,
			 const casacore::MeasurementSet& selms);

  // Local setSolve (traps lack of refant)
  virtual void setSolve(const casacore::Record& solve);
  using GJones::setSolve;

  // Default parameter value
  virtual casacore::Complex defaultPar() { return casacore::Complex(0.0); };
  
  // FIXME: was omitted
  virtual void specify(const casacore::Record& specify);
  
  // This type is not yet accumulatable
  virtual casacore::Bool accumulatable() { return false; };

  // This type is smoothable
  virtual casacore::Bool smoothable() { return true; };

  // Delay to phase calculator
  virtual void calcAllJones();

  // Hazard a guess at parameters (unneeded here)
  //  TBD?  Needed?
  virtual void guessPar(VisBuffer& ) {};

  // FringeJones uses generic gather, but solves for itself per solution
  virtual casacore::Bool useGenericGatherForSolve() { return true; };
  virtual casacore::Bool useGenericSolveOne() { return false; }

  // Post solve tinkering
  virtual void globalPostSolveTinker();

  // Local implementation of selfSolveOne (generalized signature)
  // virtual void selfSolveOne(VisBuffGroupAcc& vbga);
  virtual void selfSolveOne(SDBList&);

  virtual void solveOneVB(const VisBuffer&);

  // SolveDataBuffer is being phased out; we no longer support it.
  // virtual void solveOneSDB(const SolveDataBuffer&);

  virtual casacore::Bool& zeroRates() { return zeroRates_; }
  virtual casacore::Bool& globalSolve() { return globalSolve_; } 
  virtual casacore::Int& maxits() { return maxits_; }
  virtual casacore::Array<casacore::Double>& delayWindow() { return delayWindow_; }
  virtual casacore::Array<casacore::Double>& rateWindow() { return rateWindow_; }
  virtual casacore::Array<casacore::Bool>& paramActive() { return paramActive_; }
  virtual casacore::Bool& concatSPWs() { return concatspws_; }
  
  // Apply reference antenna
  virtual void applyRefAnt();

  virtual casacore::Int& refant() { return refant_; }
  
protected:

  // phase, delay, rate
  //  TBD:  Need to cater for parameter opt-out  (e.g., no rate solve, etc.)
  virtual casacore::Int nPar() { return 8; };

  // Jones matrix elements are NOT trivial
  virtual casacore::Bool trivialJonesElem() { return false; };

  // dJ/dp are trivial
  //  TBD: make this default in SVC?
  virtual casacore::Bool trivialDJ() { return false; };

  // Initialize trivial dJs
  //  TBD: make this default in SVC?
  virtual void initTrivDJ() {};

  // Reference frequencies
  casacore::Vector<casacore::Double> KrefFreqs_;

  
  
private:

  // Pointer to CTRateAwareTimeInterp1 factory method
  // This ensures the rates are incorporated into the time-dep interpolation
  virtual CTTIFactoryPtr cttifactoryptr() { return &CTRateAwareTimeInterp1::factory; };
  void calculateSNR(casacore::Int, casa::DelayRateFFT&);

  casacore::Int refant_; // Override
  casacore::Bool zeroRates_;
  casacore::Bool globalSolve_;
  casacore::Array<casacore::Double> delayWindow_;
  casacore::Array<casacore::Double> rateWindow_;
  casacore::Array<casacore::Bool> paramActive_;
  casacore::Int maxits_;
  casacore::Bool concatspws_;
};


} //# NAMESPACE CASA - END

#endif
