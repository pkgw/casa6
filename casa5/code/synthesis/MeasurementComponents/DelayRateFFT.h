//# DelayRateFFT.h: Declaration of fringe-fitting aux class
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

namespace casa { //# NAMESPACE CASA - BEGIN

     class SDBListGridManager {
     public:
          casacore::Double fmin_, fmax_, df_;
          casacore::Double tmin_, tmax_, dt_;
          std::set< casacore::Int > spwins_;
          std::set< casacore::Double > times_;
          casacore::Int nt_;

     protected:
          SDBList& sdbs_;
     public:
          casacore::Int nSPW() { return spwins_.size(); }
          casacore::Int getTimeIndex(casacore::Double t) { return round( (t - tmin_)/dt_ ); }
          SDBListGridManager(SDBList& sdbs) :  sdbs_(sdbs) {}



     };
     
// A utility class that provides an API that allows clients to find
// grid indices in time and frequency from a SBDList that can include
// multiple spectral windows.
     class SDBListGridManagerCombo : public SDBListGridManager {
     public:
          casacore::Int nchan_;
     public:
          std::map< casacore::Int, casacore::Int > spwPMap_; // Maps MS spws to indices in our visibility arrays
          // You can't store references in a map.
          // C++ 11 has a reference_wrapper type, but for now:
          std::map< casacore::Int, casacore::Vector<casacore::Double> const * > pspwIdToFreqMap_;
     public:
          SDBListGridManagerCombo(SDBList& sdbs_);
          // pspw is the physical spw; the one the SolveDataBuffer returns from ::spectralWindow()
          // lspw is the logical spw; the one used as an index in the fringe fitter
          casacore::Int getLSPW(casacore::Int i) { return spwPMap_.find(i)->second; }
          casacore::Int getTimeIndex(casacore::Double t) { return round( (t - tmin_)/dt_ ); }
          casacore::Float getRefFreqFromLSPW(casacore::Int lspw);
     };


     class SDBListGridManagerConcat : public SDBListGridManager {
     public:
          casacore::Double fmin_, fmax_, df_;
          casacore::Int totalChans_ = 0, nSPWChan_;
     private:
          // You can't store references in a map.
          // C++ 11 has a reference_wrapper type, but for now:
          std::map< casacore::Int, casacore::Vector<casacore::Double> const * > spwIdToFreqMap_;
     public:
          SDBListGridManagerConcat(SDBList&);
          casacore::Int bigFreqGridIndex(casacore::Double f) { return round( (f - fmin_)/df_ ); }
          casacore::Int swStartIndex(casacore::Int);
          casacore::Int nChannels() { return totalChans_;  }
          void checkAllGridpoints();
     };


// DelayRateFFT is responsible for the two-dimensional FFT of
// visibility phases to find an estimate for deay. rate and phase
// offset parameters during fringe-fitting.
     // The idiom used in KJones solvers is:
     // DelayFFT delfft1(vbga(ibuf), ptbw, refant());
     // delfft1.FFT();
     // delfft1.shift(f0[0]);
     // delfft1.searchPeak();
     // This class is designed to follow that API (without the shift).

     class DelayRateFFT {
     protected:
          casacore::Double dt_, f0_, df_;
          casacore::Int nt_;
          casacore::Int nChan_;
          casacore::Int refant_;
          // SBDListGridManager handles all the sizing and interpolating of
          // multiple spectral windows onto a single frequency grid.
          casacore::Array<casacore::Double>& delayWindow_;
          casacore::Array<casacore::Double>& rateWindow_;
          casacore::Array<casacore::Complex> Vall_;
          casacore::Array<casacore::Int> xcount_;
          casacore::Array<casacore::Float> sumw_;
          casacore::Array<casacore::Float> sumww_;
          casacore::Array<casacore::Float> peak_;
          casacore::Int nCorr_;
          casacore::Int nElem_;
          std::map< casacore::Int, std::set<casacore::Int> > activeAntennas_;
          std::set<casacore::Int> allActiveAntennas_;
          casacore::Matrix<casacore::Float> param_;
          casacore::Matrix<casacore::Bool> flag_; //?
     public:
          DelayRateFFT(SDBList& sdbs,
                       casacore::Int refant,
                       casacore::Array<casacore::Double>& delayWindow_,
                       casacore::Array<casacore::Double>& rateWindow_
               );
          virtual ~DelayRateFFT() {};
          static DelayRateFFT *makeAChild(int concat, SDBList& sdbs,
                                          casacore::Int refant,
                                          casacore::Array<casacore::Double>& delayWindow_,
                                          casacore::Array<casacore::Double>& rateWindow_);
          const casacore::Matrix<casacore::Float>& param() { return param_; }
          casacore::Int refant() const { return refant_; }
          const std::map<casacore::Int, std::set<casacore::Int> >& getActiveAntennas() const
          { return activeAntennas_; }
          const std::set<casacore::Int>& getActiveAntennasCorrelation(casacore::Int icor) const
          { return activeAntennas_.find(icor)->second; }
          void removeAntennasCorrelation(casacore::Int, std::set< casacore::Int >);
          const casacore::Matrix<casacore::Bool>& flag() const { return flag_; }
          void printActive();
          casacore::Matrix<casacore::Float> delay() const;
          casacore::Matrix<casacore::Float> rate() const;
          // If I make these NULL then the class is abstract and I can't define a factory
          // method to generate subclass instances, so I define them (recklessly!) empty
          // instead
          virtual void FFT() { 
               cerr << "DelayRateFFT::FFT()" << endl;
          } 
          virtual void searchPeak() {
               cerr << "DelayRateFFT::searchPeak()" << endl;
          }
          virtual casacore::Float snr(casacore::Int icorr, casacore::Int ielem, casacore::Float delay,
                                      casacore::Float rate) {
               cerr << "DelayRateFFTConcat::snr"<< endl;
               return 0.0; }
     };

     class DelayRateFFTCombo : public DelayRateFFT {
     private:
          SDBListGridManagerCombo gm_;
          casacore::Int nspw_;
          // 
     public:
          DelayRateFFTCombo(SDBList& sdbs, casacore::Int refant,
                            casacore::Array<casacore::Double>& delayWindow_,
                            casacore::Array<casacore::Double>& rateWindow_
               );
          DelayRateFFTCombo(casacore::Array<casacore::Complex>& data, 
                            casacore::Float f0, casacore::Float df, casacore::Float dt, SDBList& s,
                            casacore::Array<casacore::Double>& delayWindow_,
                            casacore::Array<casacore::Double>& rateWindow_
               );
          ~DelayRateFFTCombo() {};
          void FFT() override;
          void searchPeak() override;
          casacore::Float snr(casacore::Int icorr, casacore::Int ielem, casacore::Float delay,
                              casacore::Float rate) override;
          const casacore::Array<casacore::Complex>& Vall() const { return Vall_; }
          std::tuple<casacore::Double, casacore::Double, casacore::Double, casacore::Double>
               refineSearch(const casacore::Cube<casacore::Complex>&,
                            const casacore::Vector<casacore::Float>&,
                            casacore::Double, casacore::Double);
     }; // End of class DelayRateFFTCombo.

     
// DelayRateFFTConcat is responsible for the two-dimensional FFT of
// visibility phases to find an estimate for deay. rate and phase
// offset parameters during fringe-fitting.
     class DelayRateFFTConcat : public DelayRateFFT {
     private:
          SDBListGridManagerConcat gm_;
          casacore::Int nPadFactor_;
          casacore::Int nPadT_;
          casacore::Int nPadChan_;
          casacore::Double df_all_;
          casacore::Array<casacore::Complex> Vpad_;
     public:
          DelayRateFFTConcat(SDBList& sdbs, casacore::Int refant,
                             casacore::Array<casacore::Double>& delayWindow_,
                             casacore::Array<casacore::Double>& rateWindow_);
          DelayRateFFTConcat(casacore::Array<casacore::Complex>& data, casacore::Int nPadFactor,
                             casacore::Float f0, casacore::Float df, casacore::Float dt, SDBList& s,
                             casacore::Array<casacore::Double>& delayWindow_,
                             casacore::Array<casacore::Double>& rateWindow_);
          ~DelayRateFFTConcat() {};
          casacore::Double get_df_all() { return df_all_; }
          std::pair<casacore::Bool, casacore::Float> xinterp(casacore::Float alo, casacore::Float amax,
                                                              casacore::Float ahi);
          const casacore::Array<casacore::Complex>& Vpad() const { return Vpad_; }
          void FFT() override;
          void searchPeak() override;
          casacore::Float snr(casacore::Int icorr, casacore::Int ielem, casacore::Float delay,
                              casacore::Float rate) override;


     }; // End of class DelayRateFFTConcat.




} //# NAMESPACE CASA - END

