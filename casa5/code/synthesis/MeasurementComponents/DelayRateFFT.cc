//# DelayRateFFT.cc: Part of the implementation of FringeJones
//# Copyright (C) 1996,1997,1998,1999,2000,2001,2002,2003,2011
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

#include <synthesis/MeasurementComponents/FringeJones.h>

#include <msvis/MSVis/VisBuffer.h>
#include <msvis/MSVis/VisBuffAccumulator.h>
#include <ms/MeasurementSets/MSColumns.h>
#include <synthesis/CalTables/CTIter.h>
#include <synthesis/MeasurementEquations/VisEquation.h>  // *
#include <synthesis/MeasurementComponents/SolveDataBuffer.h>
#include <synthesis/MeasurementComponents/MSMetaInfoForCal.h>
#include <lattices/Lattices/ArrayLattice.h>
#include <lattices/LatticeMath/LatticeFFT.h>
#include <scimath/Mathematics/FFTServer.h>

#include <casa/Arrays/ArrayMath.h>
#include <casa/Arrays/MatrixMath.h>
#include <casa/Arrays/ArrayLogical.h>
#include <casa/BasicSL/String.h>
#include <casa/Utilities/Assert.h>
#include <casa/Exceptions/Error.h>
#include <casa/System/Aipsrc.h>

#include <casa/sstream.h>

#include <measures/Measures/MCBaseline.h>
#include <measures/Measures/MDirection.h>
#include <measures/Measures/MEpoch.h>
#include <measures/Measures/MeasTable.h>

#include <casa/Logging/LogMessage.h>
#include <casa/Logging/LogSink.h>

#include <casa/Arrays/MaskedArray.h>
#include <casa/Arrays/MaskArrMath.h>


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spblas.h>
#include <gsl/gsl_multilarge_nlinear.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_linalg.h>

#include <iomanip>                // needed for setprecision

// DEVDEBUG gates the development debugging information to standard
// error; it should be set to 0 for production.

#define DEVDEBUG false

using namespace casa::vi;
using namespace casacore;
 
namespace casa {

Complex
dotMatrixWithModel(const Matrix<Complex>& data, Double k, Double l, Double offset);

Complex
dotMatrixWithModel2(const Matrix<Complex>& f, Double k, Double l, Double offset);

Complex
dotWithOffsets(const Cube<Complex>& ft, const Vector<Float>& offsets, double k, double l);

Double
multibandFFT(const Vector<Complex>& peaks, const Vector<Float>& offsets);

tuple<Double, Double, Double, Double>
refineSearch2(const Cube<Complex>& ft,  const Vector<Float>& offsets, Double ipkt0, Double ipkch0);

tuple<Double, Double, Double, Double>
bruteForceDelay(const Cube<Complex>& ft,  const Vector<Float>& offsets, Int ipkt, Int ipkch);



static void unitize(Array<Complex>& vC) 
{
    Array<Float> vCa(amplitude(vC));
    // Divide by non-zero amps
    vCa(vCa<FLT_EPSILON)=1.0;
    vC /= vCa;
}

SDBListGridManagerCombo::SDBListGridManagerCombo(SDBList& sdbs) :
    SDBListGridManager(sdbs),
    nchan_(0)
{
    std::set<Double> fmaxes_;
    std::set<Double> fmins_;

    if (sdbs.nSDB()==0) {
        // The for loop is fine with an empty list, but code below it
        // isn't and there's nothing to lose by bailing early!
        if (DEVDEBUG) {
            cerr << "No data I guess?";
        }
        return;
    }
    
    for (Int i=0; i != sdbs_.nSDB(); i++) {
        SolveDataBuffer& sdb = sdbs_(i);
        Int pspw = sdb.spectralWindow()(0);
        Double t = sdbs_(i).time()(0);
        times_.insert(t); 
        if (spwins_.find(pspw) == spwins_.end()) {
            spwins_.insert(pspw);
            const Vector<Double>& fs = sdb.freqs();
            df_ = fs[1] - fs[0];
            pspwIdToFreqMap_[pspw] = &(sdb.freqs());
            nchan_ = max(nchan_, sdb.nChannels());
            if (DEVDEBUG) {
                cerr << "adding sdb " << i << " with " << sdb.nChannels() << " channels" << endl;
            }
            fmaxes_.insert(fs(nchan_-1));
            fmins_.insert(fs(0));
        } else {
            continue;
        }
    }
    if (nchan_ == 1) {
      nt_ = sdbs_.nSDB()/spwins_.size();
      tmin_ = *(times_.begin());
      tmax_ = *(times_.rbegin());
      dt_ = (tmax_ - tmin_)/(nt_ - 1);
      df_ = 1;
      return;
    }
    size_t i = 0;
    for (auto p=spwins_.begin(); p!=spwins_.end(); p++) {
        spwPMap_[*p] = i;
        i++;
    }
    nt_ = sdbs_.nSDB()/spwins_.size();
    tmin_ = *(times_.begin());
    tmax_ = *(times_.rbegin());
    dt_ = (tmax_ - tmin_)/(nt_ - 1);
}




Float
SDBListGridManagerCombo::getRefFreqFromLSPW(Int lspw) {
    auto p = spwins_.begin();
    std::advance(p, lspw);
    Int pspw = *p;
    const Vector<Double>& fs = *(pspwIdToFreqMap_[pspw]);
    return fs[0];
}
   
DelayRateFFT::DelayRateFFT(SDBList& sdbs, Int refant, Array<Double>& delayWindow, Array<Double>& rateWindow) :
    f0_(sdbs.centroidFreq() / 1.e9),      // GHz, for delayrate calc
    refant_(refant),
    delayWindow_(delayWindow),
    rateWindow_(rateWindow),
    Vall_(),
    xcount_(),
    sumw_(),
    sumww_(),
    activeAntennas_(),
    allActiveAntennas_()
{}

DelayRateFFT*
DelayRateFFT::makeAChild(int concat, SDBList& sdbs,
                         casacore::Int refant,
                         casacore::Array<casacore::Double>& delayWindow_,
                         casacore::Array<casacore::Double>& rateWindow_)
{
    if (concat) {
        return new DelayRateFFTConcat(sdbs, refant, delayWindow_, rateWindow_);
    } else {
        return new DelayRateFFTCombo(sdbs, refant, delayWindow_, rateWindow_);
    }
}


Matrix<Float>
DelayRateFFT::delay() const {
    IPosition start(2, 1, 0);
    IPosition stop(2, 3*nCorr_-1, nElem_-1);
    IPosition stride(2, 3, 1);
    Slicer sl(start,  stop, stride, Slicer::endIsLast);
    return param_(sl);
}

Matrix<Float>
DelayRateFFT::rate() const {
    IPosition start(2, 2, 0);
    IPosition stop(2, 3*nCorr_-1, nElem_-1);
    IPosition stride(2, 3, 1);
    Slicer sl(start,  stop, stride, Slicer::endIsLast);
    return param_(sl);
}

void DelayRateFFT::removeAntennasCorrelation(Int icor, std::set< Int > s) {
    std::set< Int > & as = activeAntennas_.find(icor)->second;
    for (std::set< Int >::iterator it=s.begin(); it!=s.end(); it++) {
        as.erase(*it);
    }
}

void
DelayRateFFT::printActive() {
    for (Int icorr=0; icorr != nCorr_; icorr++) {
        cerr << "Antennas found for correlation " << icorr << ": ";
        for (std::set<Int>::iterator it = activeAntennas_[icorr].begin();
             it != activeAntennas_[icorr].end(); it++) {
            cerr << *it << ", ";
        }
        cerr << endl;
    }
    cerr << endl;
}    



// DelayRateFFTCombo is modeled on DelayFFT in KJones.{cc|h}
DelayRateFFTCombo::DelayRateFFTCombo(SDBList& sdbs, Int refant, Array<Double>& delayWindow,
                                     Array<Double>& rateWindow) :
    DelayRateFFT(sdbs, refant, delayWindow, rateWindow),
    gm_(sdbs)
{
    nt_ = gm_.nt_;
    nChan_ = gm_.nchan_;
    nspw_ = gm_.nSPW();
    dt_ = gm_.dt_;
    df_ = gm_.df_ / 1.e9;

    if (nt_ < 2) {
        throw(AipsError("Can't do a 2-dimensional FFT on a single timestep! Please consider changing solint to avoid orphan timesteps."));
    }
    IPosition ds = delayWindow_.shape();
    if (ds.size()!=1 || ds.nelements() != 1) {
        throw AipsError("delaywindow must be a list of length 2.");
    }
    IPosition rs = rateWindow_.shape();
    if (rs.size()!=1 || rs.nelements() != 1) {
        throw AipsError("ratewindow must be a list of length 2.");
    }
    Int nCorrOrig(sdbs(0).nCorrelations());
    nCorr_ = (nCorrOrig> 1 ? 2 : 1); // number of p-hands

    for (Int i=0; i<nCorr_; i++) {
        activeAntennas_[i].insert(refant_);
    }
    
    // when we get the visCubecorrected it is already
    // reduced to parallel hands, but there isn't a
    // corresponding method for flags.
    Int corrStep = (nCorrOrig > 2 ? 3 : 1); // step for p-hands

    for (Int ibuf=0; ibuf != sdbs.nSDB(); ibuf++) {
        SolveDataBuffer& s(sdbs(ibuf));
        for (Int irow=0; irow!=s.nRows(); irow++) {
            Int a1(s.antenna1()(irow)), a2(s.antenna2()(irow));
            allActiveAntennas_.insert(a1);
            allActiveAntennas_.insert(a2);
        }
    }

    nElem_ =  1 + *(allActiveAntennas_.rbegin()) ;

    IPosition aggregateDim(2, nCorr_, nElem_, nspw_);
    xcount_.resize(aggregateDim);
    sumw_.resize(aggregateDim);
    sumww_.resize(aggregateDim);
    peak_.resize(aggregateDim);
    
    xcount_ = 0;
    sumw_ = 0.0;
    sumww_ = 0.0;
    IPosition dataSize(5, nCorr_, nElem_, nspw_, nt_, nChan_);
    Vall_.resize(dataSize);
    Int totalRows = 0;
    Int goodRows = 0;
    for (Int ibuf=0; ibuf != sdbs.nSDB(); ibuf++) {
        SolveDataBuffer& s(sdbs(ibuf));
        totalRows += s.nRows();
        if (!s.Ok())
            continue;

        Int nr = 0;
        for (Int irow=0; irow!=s.nRows(); irow++) {
            if (s.flagRow()(irow))
                continue;
            Int iant;
            Int a1(s.antenna1()(irow)), a2(s.antenna2()(irow));
            if (a1 == a2) { continue; } // we don't do autocorrelations
            else if (a1 == refant_) { iant = a2; }
            else if (a2 == refant_) { iant = a1; }
            else { continue; } // not a baseline to reference antenna
            // v has shape (nelems, ?, nrows, nchannels)
            Cube<Complex> v = s.visCubeCorrected();
            const Cube<Float>& w( s.weightSpectrum() );
            const Cube<Bool>& fl( s.flagCube() );
            Int pspw = s.spectralWindow()(0);
            Int ispw = gm_.getLSPW(pspw);
            Int t_index = gm_.getTimeIndex(s.time()(0));
            // FIXME! ispw is the Measurement Set spectral window; it
            // doesn't follow that that is addressible in our array
            // which is dimensioned for logical spectral
            // window. E.g. and i.e., if we are not joining spectral
            // windows, the spectral-window dimension of our target
            // array will always have size 1
            IPosition start (5,        0,   iant, ispw, t_index,      0);
            IPosition stop  (5,   nCorr_,      1,    1,       1, nChan_);
            IPosition stride(5,        1,      1,    1,       1,      1);
            Slicer target_slice(start, stop, stride, Slicer::endIsLength); 
            // Slicer::endIsLast is also possible
            
            Slicer source_slice(IPosition(3, 0,         0, irow),
                                IPosition(3, nCorr_,  nChan_, 1),
                                IPosition(3, corrStep,      1,  1), Slicer::endIsLength);
                
            Slicer flagSlice(IPosition(3, 0,             0, irow),
                             IPosition(3, nCorr_,    nChan_, 1),
                             IPosition(3, corrStep,       1, 1), Slicer::endIsLength);
            nr++;
            Array<Complex>rhs( v(source_slice).nonDegenerate(1) );
            const Array<Float>& weights( w(source_slice).nonDegenerate(1) );
            unitize(rhs); 
            Vall_(target_slice).nonDegenerate(1) = rhs * weights;
            const Array<Bool>& flagged( fl(flagSlice).nonDegenerate(1) );

            // Zero flagged entries.
            Vall_(target_slice).nonDegenerate(1)(flagged) = Complex(0.0);

            if (!allTrue(flagged)) {
                for (Int icorr=0; icorr<nCorr_; ++icorr) {
                    IPosition p(2, icorr, iant);
                    Bool actually = false;
                    activeAntennas_[icorr].insert(iant);
                    for (Int ichan=0; ichan != (Int) nChan_; ichan++) {
                        IPosition pchan(2, icorr, ichan);
                        if (!flagged(pchan)) {
                            Float wv = weights(pchan);
                            xcount_(p)++;
                            sumw_(p) += wv;
                            sumww_(p) += wv*wv;
                            actually = true;
                        }
                    }
                    if (actually) {
                        activeAntennas_[icorr].insert(iant);
                        goodRows++;
                    }
                }
            }
        }
    }
}


// What even *is* this?
DelayRateFFTCombo::DelayRateFFTCombo(Array<Complex>& data, Float f0, Float df, Float dt, SDBList& s,
                                     Array<Double>& delayWindow, Array<Double>& rateWindow) :
    DelayRateFFT(s, -1, delayWindow, rateWindow),
    gm_(s)
{
    dt_ = dt;
    df_ = df;
    IPosition shape = data.shape();
    nCorr_ = shape(0);
    nElem_ = shape(1);
    nspw_ = shape(2);
    nt_ = shape(3);
    nChan_ = shape(4);
    IPosition dataSize(5, nCorr_, nElem_, nspw_, nt_, nChan_);
    Vall_.resize(dataSize);
    
    IPosition start(5, 0, 0, 0, 0, 0);
    IPosition stop(5, nCorr_,  nElem_, nspw_, nt_, nChan_);
    IPosition stride(5, 1, 1, 1, 1);
    Slicer target_slice(start, stop, stride, Slicer::endIsLength);
    Vall_(target_slice) = data;

    unitize(Vall_);

}

void
DelayRateFFTCombo::FFT() {
    if (DEVDEBUG) {
        cerr << "DelayRateFFTCombo::FFT()" << endl;
    }
    // Axes are 0: correlation (i.e., hand of polarization), 1: antenna, 2: spectral window, 3: time, 4: channel
    // the machinery is there to say that we only want to FFT the last two axes
    Vector<Bool> ax(5, false);
    ax(3) = true;
    ax(4) = true;
    // Also copied from DelayFFT in KJones.
    // we make a copy to FFT in place.
    ArrayLattice<Complex> c(Vall_);
    LatticeFFT::cfft0(c, ax, true);
    if (DEVDEBUG) {
        cerr << "FFT transformed" << endl;
    }
}

// In the new paradigm, this is where the new stuff happens. Instead of
// interpolating the peaks on a single big grid, we have to synthesise
// our estimate from separately FFT-ed spectral windows, using the new
// off-grid peak formalism
void
DelayRateFFTCombo::searchPeak() {
    if (DEVDEBUG) {
        cerr << "DelayRateFFTCombo::searchPeak()" << endl;
    }
    
    // Recall param_ -> [phase, delay, rate] for each correlation
    param_.resize(3*nCorr_, nElem_); // This might be better done elsewhere.
    param_.set(0.0);
    flag_.resize(3*nCorr_, nElem_);
    flag_.set(true);  // all flagged initially

    Double bw = Float(nChan_)*df_;
    for (Int icorr=0; icorr<nCorr_; ++icorr) {
        flag_(icorr*3 + 0, refant()) = false; 
        flag_(icorr*3 + 1, refant()) = false;
        flag_(icorr*3 + 2, refant()) = false;
        for (Int ielem=0; ielem<nElem_; ++ielem) {
            if (ielem==refant()) {
                continue;
            }
            if (activeAntennas_[icorr].find(ielem)==activeAntennas_[icorr].end()) {
                continue;
            }

            // Below is the gory details for turning delay window into index range
            Int sgn = (ielem < refant()) ? 1 : -1;
            Double d0 = sgn*delayWindow_(IPosition(1, 0));
            Double d1 = sgn*delayWindow_(IPosition(1, 1));
            if (d0 > d1) std::swap(d0, d1);
            d0 = max(d0, -0.5/df_);
            d1 = min(d1, (0.5-1/nChan_)/df_);
    
            // It's simpler to keep the ranges as signed integers and
            // handle the wrapping of the FFT in the loop over
            // indices. Recall that the FFT result returned has indices
            // that run from 0 to nChan_/2 -1 and then from
            // -nChan/2 to -1, so far as our delay is concerned.
            Int i0 = bw*d0;
            Int i1 = bw*d1;
            if (i1==i0) i1++;
            // Now for the gory details for turning rate window into index range
            Double width = nt_*dt_*1e9*f0_;
            Double r0 = sgn*rateWindow_(IPosition(1,0));
            Double r1 = sgn*rateWindow_(IPosition(1,1));
            if (r0 > r1) std::swap(r0, r1);
            r0 = max(r0, -0.5/(dt_*1e9*f0_));
            r1 = min(r1, (0.5 - 1/nt_)/(dt_*1e9*f0_));
    
            Int j0 = width*r0;
            Int j1 = width*r1;
            if (j1==j0) j1++;
            // FIXME: We now want an incoherent sum of the amplitudes of all the subbands!
            Matrix<Float> inco(IPosition(2, nt_, nChan_));
            inco = 0; // Note, painfully, that this is not the default!
            // NB: Time, Channel
            // And once again we fail at slicing
            // IPosition stop (5,     1,     1,  spw_, nt_, nChan_);
            // FIXME: we shouldn't but we will: just choose the first spw
            for (size_t ispw=0; ispw!=size_t(nspw_); ispw++) {
                IPosition start(5, icorr, ielem,     ispw,   0,      0);
                IPosition stop (5,     1,     1,     1,    nt_, nChan_);
                IPosition step (5,     1,     1,     1,      1,      1);
                Slicer sl(
                    start,
                    stop,
                    step,
                    Slicer::endIsLength);
                const Matrix<Complex>& aS(Vall_(sl).nonDegenerate());
                inco += amplitude(aS);
            }

            // We search the incoherent sum for the position of the
            // maximum on the grid, which we will use to start the
            // refinement process
            Int ipkch(0);
            Int ipkt(0);
            Float amax(-1.0);
            for (Int itime0=j0; itime0 != j1; itime0++) {
                Int itime = (itime0 < 0) ? itime0 + nt_ : itime0;
                for (Int ich0=i0; ich0 != i1; ich0++) {
                    Int ich = (ich0 < 0) ? ich0 + nChan_ : ich0;
                    if (inco(itime, ich) > amax) {
                        ipkch = ich;
                        ipkt  = itime;
                        amax=inco(itime, ich);
                    }
                }
            }

            Float phase0;
            for (size_t ispw=0; ispw!=size_t(nspw_); ispw++) {
                Complex p = Vall_(IPosition(5, icorr, ielem, ispw, ipkt, ipkch));
                if (ispw==0) {
                    phase0 = arg(p);
                }
                if (DEVDEBUG) {
                    cerr << "   Before refining: " << "ispw " << ispw << " peak " << abs(p) << " ang " << arg(p) << endl;
                }
            }
            
            // Finished grovelling. Now we have the location of the
            // maximum amplitude on the grid, and we refine it by
            // looking off the grid a little
            Array<Complex> blVis(
                Vall_(Slicer(
                          IPosition(5, icorr, ielem,     0,   0,      0),
                          IPosition(5,     1,     1, nspw_, nt_, nChan_),
                          IPosition(5,     1,     1,     1,   1,      1),
                          Slicer::endIsLength)).nonDegenerate(IPosition(1,2)));

            Vector<Float> offsets(nspw_);
            Double bw = gm_.nchan_ * gm_.df_;
            for (size_t lspw=0; lspw!=size_t(nspw_); lspw++) {
                Double f = gm_.getRefFreqFromLSPW(lspw);
                Double f0 = gm_.getRefFreqFromLSPW(0);
                offsets(lspw) = (f - f0)/bw;
            }
            if (DEVDEBUG) {
                cerr << "Offsets: " << offsets << endl;
            }
            Complex c0 = dotWithOffsets(blVis, offsets, double(ipkt), double(ipkch));

            Vector<Complex> peaks(
                blVis(Slicer(
                          IPosition(3, 0,  ipkt, ipkch),
                          IPosition(3, nspw_, 1, 1),
                          IPosition(3, 1, 1, 1),
                          Slicer::endIsLength)).nonDegenerate(1));
            
            // In units of spw BW, like offsets
            Double dpkch = multibandFFT(peaks, offsets);
            if (DEVDEBUG) {
                cerr << "Peaks: " << peaks << endl;
                cerr << "dpkch " << dpkch << endl;
            }
            // FIXME: Sign convention dilemma: do we add or subtract dpkch?
            tuple<Double, Double, Double, Double> p = refineSearch(blVis, offsets, ipkt, ipkch + dpkch);
            Double pkt   = std::get<0>(p);
            Double pkch  = std::get<1>(p);
            Double peak  = std::get<2>(p);
            Double phase = std::get<3>(p);
            if (DEVDEBUG) {
                cerr << "[DelayRateFFTCombo::SearchPeak] " << "icorr " << icorr << " ielem " << ielem
                     << " from (" << ipkt << ", " << ipkch << ", peak "  << abs(c0) << " (incoherently " << inco(ipkt, ipkch) << "), angle " << phase0 << ")" 
                     << " to (" << pkt << ", " << pkch << ", peak " << peak << ", angle " << phase <<  ")"
                     << endl;
            }
            peak_(IPosition(2, icorr, ielem)) = peak;
            param_(icorr*3 + 0, ielem) = sgn*phase;
            Float delay = (pkch)/Float(nChan_);
            if (delay > 0.5) delay -= 1.0;           // fold
            delay /= df_;                            // nsec
            param_(icorr*3 + 1, ielem) = sgn*delay; 
            Double rate = (pkt)/Float(nt_);
            if (rate > 0.5) rate -= 1.0;
            Double rate0 = rate/dt_;
            Double rate1 = rate0/(1e9 * f0_); 
            param_(icorr*3 + 2, ielem) = Float(sgn*rate1);
            if (DEVDEBUG) {
                cerr << "delay " << sgn*delay << " rate1 " << sgn*rate1 << endl;
            }
            // Set 3 flags.
            flag_(icorr*3 + 0, ielem)=false; 
            flag_(icorr*3 + 1, ielem)=false;
            flag_(icorr*3 + 2, ielem)=false;
            //cerr << "Set everything " << endl;
        }
    }
}

Double
multibandFFT(const Vector<Complex>& peaks, const Vector<Float>& offsets) {
    // We take the individual band phases from the peaks of the SPWs
    // and assume they are the peculiar phases for their SPW; then we
    // calculate the Direct Discrete FT of those phases on their
    // offsets and read the peak off of that
    Int nspw = peaks.size();

    Vector<Float> amps( amplitude(peaks) );
    Int len;
    amps.shape(len);

    Float eps = 1e-10;
    std::set<size_t> badInds;
    for (size_t i=0; i!=size_t(len); i++) {
        if (amps[i] < eps) {
            if (DEVDEBUG) {
                cerr << " blacklisting " << i << endl;
            }
            badInds.insert(i);
        }
    }
    Vector<Complex> phasors( peaks/amplitude(peaks) );

    // I want to print angles but I can't get the API to do that
    // in a reasonable way; arg(phasors) doesn't work
    // Vector<Float> angles ( arg(phasors) );
    if (DEVDEBUG) {
        cerr << "phasors " << phasors << endl;
    }
    // FIXME: the last offset is at the *beginning* of the subband!
    // size_t nbins=2*size_t(max(offsets)+1+0.5);
    size_t nbins=size_t(max(offsets)+1+0.5);

    Vector<Float> freqs(nbins);
    freqs = 0;
    for (size_t k=0; k!=nbins; k++) {
        freqs[k] = (Float(k) - nbins/2)/float(nbins);
    }
    Vector<Complex> X(nbins);
    X = 0;
    for (size_t n=0; n!=size_t(nspw); n++) {
        if (badInds.find(n) != badInds.end()) {
            if (DEVDEBUG) {
                cerr << " skipping " << n << endl;
            }
            continue;
        }
        for (size_t k=0; k!=nbins; k++) {
            // FIXME: also swap signs here!
            X[k] += phasors[n]*exp(-Complex(0,1)*Complex(C::_2pi*offsets(n)*freqs(k)));
        }
    }
    Int k_max = -1;
    Float zmax = -1.0;
    for (size_t k=0; k!=nbins; k++) {
        Complex z = X[k];
        Float a = abs(z);
        if (DEVDEBUG) {
            cerr << "k " << k << " freq " << freqs(k) << " a " << a << endl;
        }
        if (a>zmax) {
            k_max = k;
            zmax = a;
        }
    }
    if (DEVDEBUG) {
        cerr << "k_max = " << k_max << endl;
    }
    return freqs(k_max);
}

// Auxilliary function
Complex
dotWithOffsets(const Cube<Complex>& ft, const Vector<Float>& offsets, double k, double l) {
    Int nspw = ft.nrow();
    Int ni = ft.ncolumn();
    Int nj = ft.nplane();

    // cerr << "ni " << ni << " nj " << nj << endl;
    if (k<0) k += (ni);
    if (l<0) l += (nj);
    Complex p(0.0, 0.0);
    for (Int s=0; s!=nspw; s++) {
        const Matrix<Complex>& ft_s = ft.yzPlane(s);
        Double f_off = offsets(s);
        Complex r = dotMatrixWithModel(ft_s, k, l, f_off);
        p += r;
    }    
    return p;
}

// We need a function to minimize, which has to return a real value,
// but when we're done we need to find the peak and also its argument,
// so we factor out the complex version...
Complex
c_peak_fn(const gsl_vector *x, void *vparams) {
    std::pair< Cube<Complex> const *, Vector<Float> const * > *params =
        (std::pair< Cube<Complex> const *, Vector<Float> const * > *)(vparams );

    double k = gsl_vector_get(x, 0);
    double l = gsl_vector_get(x, 1);

    const Cube<Complex>& ft ( *(params->first ));
    const Vector<Float>& offsets ( *(params->second));
    return dotWithOffsets(ft, offsets, k, l);
}

// ... and then call it from a real version that we can give to the
// optimizer (with a sign change; by tradition these are minimizer
// routines and we want a maximum)
double
my_peak_fn(const gsl_vector *x, void *vparams) {
    Complex p = c_peak_fn(x, vparams);
    return -abs(p);
}

tuple<Double, Double, Double, Double>
bruteForceDelay(const Cube<Complex>& ft,  const Vector<Float>& offsets, Int ipkt, Int ipkch) {
    Double k(ipkt);
    Double l(ipkch);
    Int n = 10;

    Double l_p = l;
    Complex peak = Complex(0, 0);
    for (Int i=-n; i!=n+1; i++) {
        Double dch = double(i)/(2*n);
        Complex p = dotWithOffsets(ft, offsets, k, l+dch);
        if (abs(p) > abs(peak)) {
            if (DEVDEBUG) {
                cerr << "[bruteForceDelay] l " << l+dch << " peak " << abs(p) << endl;
            }
            peak = p;
            l_p = l+dch;
            
        }
    }            
    tuple<Double, Double, Double, Double> t;
    t = std::make_tuple(k, l_p, abs(peak), arg(peak));
    return t;
}


tuple<Double, Double, Double, Double>
DelayRateFFTCombo::refineSearch(const Cube<Complex>& ft,  const Vector<Float>& offsets, Double pkt, Double pkch) {
    // small@jive.eu: I borrowed most of this code from the GSL documentation of multimin:
    // <https://www.gnu.org/software/gsl/doc/html/multimin.html>
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
    /* Starting point */
    gsl_vector *x = gsl_vector_alloc (2);
    gsl_vector_set(x, 0, pkt);
    gsl_vector_set(x, 1, pkch);
    if (DEVDEBUG) {
        cerr << "pkch starts at " << pkch << endl;
    }
    /* Set initial step sizes */
    /* Fixme! I need to think harder about this value! */
    gsl_vector* steps = gsl_vector_alloc (2);
    gsl_vector_set_all(steps, 0.005);
    
    /* Initialize method and iterate */
    
    gsl_multimin_function minex_func;
    std::pair< Cube<Complex> const * , Vector<Float> const * > par = make_pair(&ft, &offsets);
     
    minex_func.n = 2;
    minex_func.f = my_peak_fn;
    minex_func.params = &par;
     
    gsl_multimin_fminimizer *s = gsl_multimin_fminimizer_alloc(T, 2);
    gsl_multimin_fminimizer_set (s, &minex_func, x, steps);


    int status;
    size_t iter = 0;
    do {
        iter++;
        status = gsl_multimin_fminimizer_iterate(s);
        if (status) break;
        double size = gsl_multimin_fminimizer_size(s);
        status = gsl_multimin_test_size(size, 1e-4);
        if (status == GSL_SUCCESS) {
            printf ("converged to minimum at\n");
        }
        printf("%5d ipkt %10.3e ipkch %10.3e f() = %7.3f size = %10.3f\n",
               iter,
               gsl_vector_get(s->x, 0),
               gsl_vector_get(s->x, 1),
               s->fval,
               size);
    } while (status == GSL_CONTINUE && iter < 2);
    // while (status == GSL_CONTINUE && iter < 100);
    tuple<Double, Double, Double, Double> p;
    // if (status == GSL_SUCCESS) {
    if (1) {
        // Double peak = gsl_multimin_fminimizer_minimum(s);
        Double ipkt  = gsl_vector_get(s->x, 0);
        Double ipkch = gsl_vector_get(s->x, 1);
        Complex c = c_peak_fn(s->x, &par);
        p = std::make_tuple(ipkt, ipkch, abs(c), arg(c));
    } else {
        // FIXME: More spurious zeros!
        p = std::make_tuple(pkt, pkch, 0.0, 0.0);
    }
    gsl_vector_free(x);
    gsl_vector_free(steps);
    gsl_multimin_fminimizer_free(s);
     
    return p;
}

Float
DelayRateFFTCombo::snr(Int icorr, Int ielem, Float delay, Float rate) {
    if (DEVDEBUG) {
        cerr << "DelayRateFFTCombo::snr"<< endl;
    }
    // We calculate a signal-to-noise ration for the 2D FFT fringefit
    // using a formula transcribed from AIPS FRING.
    //
    // Have to convert delay and rate back into indices on the padded 2D grid.
    IPosition p(2, icorr, ielem);
    Float peak = peak_(p);
    if (peak > 0.999*sumw_(p)) {
        cerr << "Clipping peak for element " << ielem << " from " << peak << " to " << 0.999*sumw_(p) << endl;
        peak=0.999*sumw_(p);
    }
    // xcount is number of data points for baseline to ielem
    // sumw is sum of weights,
    // sumww is sum of squares of weights
    Float cwt;
    if (fabs(sumw_(p))<FLT_EPSILON) {
        cwt = 0;
        if (DEVDEBUG) {
            cerr << "Correlation " << icorr << " antenna " << ielem << ": sum of weights is zero." << endl;
        }
    } else {
        Float x = C::pi/2*peak/sumw_(p);
        // The magic numbers in the following formula are from AIPS FRING
        cwt = (pow(tan(x), 1.163) * sqrt(sumw_(p)/sqrt(sumww_(p)/xcount_(p))));
        if (DEVDEBUG) {
            cerr << "Correlation " << icorr << " antenna " << ielem 
                 << " peak=" << peak << "; xang=" << x << "; xcount=" << xcount_(p) << "; sumw=" << sumw_(p) << "; sumww=" << sumww_(p)
                 << " snr " << cwt << endl;
        }
    }
    return cwt;
}
    

// There isn't a usable sinc lying around that I can see
Double sinc(Double x) {
    Double p = C::pi*x;
    return sin(p)/p;
}

// This implements a Fancy Sinc strategy for finding the peak of a
// single SPW. I no longer know if this can be extended to multiple
// SPWs.
tuple<Double, Double, Double, Double>
refineSearch2(const Cube<Complex>& ft, const Vector<Float>& offsets, Int ipkt0, Int ipkch0) {
    Double imax, jmax;
    // FIXME! We need to do a single SPW first
    if (ft.nrow() > 1) {
        throw AipsError("Only one spw for now");
    }
    Int s = 0;
    const Matrix<Complex>& ft_s(ft.yzPlane(s));
    size_t ni = ft_s.nrow();
    size_t nj = ft_s.ncolumn();
    // we need to do a roll left/down and sum but Casacore doesn't have a matrix roll so we may as well do it by hand
    Float pmax = -1;
    size_t ipkt(0);
    size_t ipkch(0);
    Matrix<Float> sumAbs(ni, nj);
    for (size_t i=0; i!=ni; i++) {
        // FIXME: The variable i1 is unused!
        // Either use it or lose it!
        // size_t i1 = (i==ni-1) ? 0 : i+1;
        for (size_t j=0; j!=nj; j++) {
            size_t j1 = (j==nj-1) ? 0 : j+1;
            Float sumAbs = (abs(ft_s(i, j)) + abs(ft_s(i, j1)));
            if ((sumAbs) > pmax) {
                ipkt = i;
                ipkch = j;
                pmax = sumAbs;
            }
        }
    }
    Double y0 = abs(ft_s(ipkt, ipkch));
    Double y1 = abs(ft_s(ipkt, (ipkch==nj-1)? 0 : ipkch+1));
    Double ym1 = abs(ft_s(ipkt, (ipkch==0)? nj-1 : ipkch-1));

    Double z1 =  abs(ft_s((ipkt==ni-1) ? 0 : ipkt+1, ipkch));
    Double zm1 =  abs(ft_s((ipkt==0) ? ni-1 : ipkt-1, ipkch));
    
    
    Double d0 = y1/(y0+y1);
    Double peak( y0/sinc(d0) );
    if (DEVDEBUG) {
        cerr << "y0 " << y0 << " y1 " << y1 << " ym1 " << ym1 << " d0 " << d0 << " sinc(d0) "
             << sinc(d0) << " peak " << peak << endl;
        cerr << "y0 " << y0 << " z1 " << z1 << " zm1 " << zm1 << endl;
        cerr << real(abs(ft_s.row(ipkt)))/y0 << endl;
    }
    Double phase( 0.0 );
    imax = Double(ipkt);
    jmax = Double(ipkch) + d0;
    tuple<Double, Double, Double, Double> p = std::make_tuple(imax, jmax, peak, phase);
    return p;
}

Complex
dotMatrixWithModel(const Matrix<Complex>& data, Double k, Double l, Double offset)
{
    size_t ni = data.nrow();
    size_t nj = data.ncolumn();
    Matrix<Complex> model(ni, nj);

    // Note that k is the time index of the array, l is the frequency index.
    // Offsets correspond to frequencies
    Double eps = 1e-8;
    Int k_int = floor(k);
    Double k_del = k - k_int;
    Bool k_flag =  (fabs(k_del) < eps);
    Int l_int = floor(l);
    Double l_del = l - l_int;
    Bool l_flag = (fabs(l_del) < eps);
    Complex t0;
    Complex t1;
    // FIXME! We're messing with reversing signs
    for (size_t i=0; i!=ni; i++) {
        if (k_flag) {
            t0 = Complex(i==k);
        } else { // if k isn't an integer!
            t0 = ( (1-exp(Complex(0, -1.0*C::_2pi*(k-i))))/
                   (1-exp(Complex(0, -1.0*C::_2pi*(k-i)/Double(ni)))) );
            t0 /= ni;
        }
        for (size_t j=0; j!=nj; j++) {
            if (l_flag) {
                t1 = Complex(j==l);
            } else { // if l isn't an integer!
                t1 = ( (1-exp(Complex(0, -1.0*C::_2pi*(l-j))))/
                       (1-exp(Complex(0, -1.0*C::_2pi*(l-j)/Double(nj)))) );
                t1 /= nj;
            }
            model(i, j) = exp(Complex(0, -1.0*C::_2pi*offset*l_del))*t0*t1;
        }
    }
    Complex t2 = sum(data*model);
    return t2;
}


Complex
dotMatrixWithModel2(const Matrix<Complex>& data, Double k0, Double l0, Double offset)
{
    size_t ni = data.nrow();
    size_t nj = data.ncolumn();
    Matrix<Complex> model(ni, nj);
    model = 0.0;
    
    Double eps = 1e-8;
    Int k_int = Int(ceil(k0));
    Int l_int = Int(ceil(l0));
    
    Double d0 = k_int - k0;
    Double d1 = l_int - l0;

    Bool k_flag = (fabs(d0) < eps);
    Bool l_flag = (fabs(d1) < eps);

    Complex c0, c1;

    Complex s (0.0);
    Complex t0 = Complex(sin(C::pi*d0)/C::pi)*exp(Complex(0, C::pi*d0));
    Complex t1 = Complex(sin(C::pi*d1)/C::pi)*exp(Complex(0, C::pi*d1));
    for (Int dk=-2; dk!=2; dk++) {
        size_t k = (k_int + dk + ni) % ni;
        if (k_flag) {
            c0 = Complex(dk==0);
        } else { // if k isn't an integer!
            c0 = Complex(1/(d0+dk))*t0;
        }
        for (Int dl=-2; dl!=2; dl++) {
            size_t l = (l_int+dl+nj) % nj;
            if (l_flag) {
                c1 = Complex(dl==0);
            } else { // if l isn't an integer!
                c1 = Complex(1/(d1+dl))*t1;
            }
            Complex d = data(k, l);
            Complex m = c0*c1;
            s += d*m;
        }
    }
    return s;
}


SDBListGridManagerConcat::SDBListGridManagerConcat(SDBList& sdbs) :
    SDBListGridManager(sdbs)
{
    std::set<Double> fmaxes;
    std::set<Double> fmins;
    Float dfn(0.0);
    Int totalChans0(0) ;
    Int nchan(0);

    if (sdbs_.nSDB()==0) {
        // The for loop is fine with an empty list, but code below it
        // isn't and there's nothing to lose by bailing early!
        return;
    }
        
    for (Int i=0; i != sdbs_.nSDB(); i++) {
        SolveDataBuffer& sdb = sdbs_(i);
        Int spw = sdb.spectralWindow()(0);
        Double t = sdbs_(i).time()(0);
        times_.insert(t); 
        if (spwins_.find(spw) == spwins_.end()) {
            spwins_.insert(spw);
            const Vector<Double>& fs = sdb.freqs();
            spwIdToFreqMap_[spw] = &(sdb.freqs());
            nchan = sdb.nChannels();
            fmaxes.insert(fs(nchan-1));
            fmins.insert(fs(0));
            // We assume they're all at the same time.

            totalChans0 += nchan;
            Float df0 = fs(1) - fs(0);
            dfn = (fs(nchan-1) - fs(0))/(nchan-1);
            if (DEVDEBUG) {
                cerr << "Spectral window " << spw << " has " << nchan << " channels" << endl;
                cerr << "df0 "<< df0 << "; " << "dfn " << dfn << endl;
            }
        } else {
            continue;
        }
    }
    if (nchan == 1) {
      nt_ = sdbs_.nSDB()/spwins_.size();
      tmin_ = *(times_.begin());
      tmax_ = *(times_.rbegin());
      dt_ = (tmax_ - tmin_)/(nt_ - 1);
      nSPWChan_ = nchan;
      fmin_ = *(fmins.begin());
      fmax_ = fmin_;
      totalChans_ = 1;
      df_ = 1;
      return;
    }
      
    nt_ = sdbs_.nSDB()/spwins_.size();
    tmin_ = *(times_.begin());
    tmax_ = *(times_.rbegin());
    dt_ = (tmax_ - tmin_)/(nt_ - 1);
    nSPWChan_ = nchan;
    fmin_ = *(fmins.begin());
    fmax_ = *(fmaxes.rbegin());
    totalChans_ = round((fmax_ - fmin_)/dfn + 1);
    df_ = (fmax_ - fmin_)/(totalChans_-1);
    if (DEVDEBUG) {
        cerr << "Global fmin " << fmin_ << " global max " << fmax_ << endl;
        cerr << "nt " << nt_ << " dt " << dt_ << endl;
        cerr << "tmin " << tmin_ << " tmax " << tmax_ << endl;
        cerr << "Global dt " << tmax_ - tmin_ << endl;
        cerr << "Global df " << df_ << endl;
        cerr << "I guess we'll need " << totalChans_ << " freq points in total." << endl;
        cerr << "Compared to " << totalChans0 << " with simple-minded concatenation." << endl;
        cerr << "spwins_.size() " << spwins_.size() << endl;
        cerr << "nSPW() " << nSPW() << endl;
    }
}

// checkAllGridPoints is a diagnostic funtion that should not be called
// in production releases, but it doesn't do any harm to have it
// latent.
void
SDBListGridManagerConcat::checkAllGridpoints() {
    map<Int, Vector<Double> const *>::iterator it;
    for (it = spwIdToFreqMap_.begin(); it != spwIdToFreqMap_.end(); it++) {
        Int spwid = it->first;
        Vector<Double> const* fs = it->second;
        Int length;
        fs->shape(length);
        for (Int i=0; i!=length; i++) {
            Double f = (*fs)(i);
            Int j = bigFreqGridIndex(f);
            if (DEVDEBUG) {
                cerr << "spwid, i = (" << spwid << ", " << i << ") => " << j << " (" << f << ")" << endl;
            }
        }
    }
    if (DEVDEBUG) {
        cerr << "[1] spwins.size() " << nSPW() << endl;
    }
}

Int
SDBListGridManagerConcat::swStartIndex(Int spw) {
    Vector<Double> const* fs = spwIdToFreqMap_[spw];
    Double f0 = (*fs)(0);
    return bigFreqGridIndex(f0);
}
   
DelayRateFFTConcat::DelayRateFFTConcat(SDBList& sdbs, Int refant, Array<Double>& delayWindow,
                                       Array<Double>& rateWindow) :
    DelayRateFFT(sdbs, refant, delayWindow, rateWindow),
    gm_(sdbs)
{
    if (DEVDEBUG) {
        cerr << "gm_.nSPW() " << gm_.nSPW() << endl;
    }
    nPadFactor_= max(2, 8  / gm_.nSPW());
    nt_ = gm_.nt_;
    nPadT_ = nPadFactor_ * nt_;
    nChan_ = gm_.nChannels();
    nPadChan_ = nPadFactor_*nChan_;
    dt_ = gm_.dt_;
    f0_ = sdbs.centroidFreq() / 1.e9;      // GHz, for delayrate calc
    df_ = gm_.df_ / 1.e9;
    df_all_ = gm_.fmax_ - gm_.fmin_;
    // This check should be commented out in production:
    // gm_.checkAllGridpoints();
    if (nt_ < 2) {
        throw(AipsError("Can't do a 2-dimensional FFT on a single timestep! Please consider changing solint to avoid orphan timesteps."));
    }
    IPosition ds = delayWindow_.shape();
    if (ds.size()!=1 || ds.nelements() != 1) {
        throw AipsError("delaywindow must be a list of length 2.");
    }
    IPosition rs = rateWindow_.shape();
    if (rs.size()!=1 || rs.nelements() != 1) {
        throw AipsError("ratewindow must be a list of length 2.");
    }
    Int nCorrOrig(sdbs(0).nCorrelations());
    nCorr_ = (nCorrOrig> 1 ? 2 : 1); // number of p-hands

    for (Int i=0; i<nCorr_; i++) {
        activeAntennas_[i].insert(refant_);
    }
    
    // when we get the visCubecorrected it is already
    // reduced to parallel hands, but there isn't a
    // corresponding method for flags.
    Int corrStep = (nCorrOrig > 2 ? 3 : 1); // step for p-hands
    for (Int ibuf=0; ibuf != sdbs.nSDB(); ibuf++) {
        SolveDataBuffer& s(sdbs(ibuf));
        for (Int irow=0; irow!=s.nRows(); irow++) {
            Int a1(s.antenna1()(irow)), a2(s.antenna2()(irow));
            allActiveAntennas_.insert(a1);
            allActiveAntennas_.insert(a2);
        }
    }
    nElem_ =  1 + *(allActiveAntennas_.rbegin()) ;
    IPosition aggregateDim(2, nCorr_, nElem_);
    xcount_.resize(aggregateDim);
    sumw_.resize(aggregateDim);
    sumww_.resize(aggregateDim);

    xcount_ = 0;
    sumw_ = 0.0;
    sumww_ = 0.0;
    
    if (DEVDEBUG) {
        cerr << "Filling FFT grid with " << sdbs.nSDB() << " data buffers." << endl;
    }
    IPosition paddedDataSize(4, nCorr_, nElem_, nPadT_, nPadChan_);
    Vpad_.resize(paddedDataSize);
    Int totalRows = 0;
    Int goodRows = 0;

    
    for (Int ibuf=0; ibuf != sdbs.nSDB(); ibuf++) {
        SolveDataBuffer& s(sdbs(ibuf));
        totalRows += s.nRows();
        if (!s.Ok())
            continue;

        Int nr = 0;
        for (Int irow=0; irow!=s.nRows(); irow++) {
            if (s.flagRow()(irow))
                continue;
            Int iant;
            Int a1(s.antenna1()(irow)), a2(s.antenna2()(irow));
            if (a1 == a2) {
                continue;
            }
            else if (a1 == refant_) {
                iant = a2;
            }
            else if (a2 == refant_) {
                iant = a1;
            }
            else {
                continue;
            }
            // OK, we're not skipping this one so we have to do something.

            // v has shape (nelems, ?, nrows, nchannels)
            Cube<Complex> v = s.visCubeCorrected();
            Cube<Float> w = s.weightSpectrum();
            Cube<Bool> fl = s.flagCube();
            Int spw = s.spectralWindow()(0);
            Int f_index = gm_.swStartIndex(spw);    // ditto!
            Int t_index = gm_.getTimeIndex(s.time()(0));
            Int spwchans = gm_.nSPWChan_;
            IPosition start(4,      0,      iant, t_index, f_index);
            IPosition stop(4,      nCorr_,    1,       1, spwchans);
            IPosition stride(4,      1,         1,       1, 1);
            Slicer sl1(start,     stop, stride, Slicer::endIsLength);
            Slicer sl2(IPosition(3, 0,         0, irow),
                       IPosition(3, nCorr_, spwchans, 1),
                       IPosition(3, corrStep,        1,  1), Slicer::endIsLength);
                
            Slicer flagSlice(IPosition(3, 0,         0, irow),
                             IPosition(3, nCorr_, spwchans, 1),
                             IPosition(3, corrStep,        1, 1), Slicer::endIsLength);
            nr++;
            if (DEVDEBUG) {
                cerr << "nr " << nr
                     << " irow " << endl
                     << "Vpad shape " << Vpad_.shape() << endl
                     << "v shape " << v.shape() << endl
                     << "sl2 " << sl2 << endl
                     << "sl1 " << sl1 << endl
                     << "flagSlice " << flagSlice << endl;
            }
            Array<Complex> rhs = v(sl2).nonDegenerate(1);
            Array<Float> weights = w(sl2).nonDegenerate(1);
                
            unitize(rhs);
            Vpad_(sl1).nonDegenerate(1) = rhs * weights;

            Array<Bool> flagged(fl(flagSlice).nonDegenerate(1));
            // Zero flagged entries.
            Vpad_(sl1).nonDegenerate(1)(flagged) = Complex(0.0);

            if (!allTrue(flagged)) {
                for (Int icorr=0; icorr<nCorr_; ++icorr) {
                    IPosition p(2, icorr, iant);
                    Bool actually = false;
                    activeAntennas_[icorr].insert(iant);
                    for (Int ichan=0; ichan != (Int) spwchans; ichan++) {
                        IPosition pchan(2, icorr, ichan);
                        if (!flagged(pchan)) {
                            Float wv = weights(pchan);
                            if (wv < 0) {
                                cerr << "spwchans " << spwchans << endl;
                                cerr << "Negative weight << (" << wv << ") on row "
                                     << irow << " baseline (" << a1 << ", " << a2 << ") "
                                     << " channel " << ichan << endl;
                                cerr << "pchan " << pchan << endl;
                                cerr << "Weights " << weights << endl;
                            }
                            xcount_(p)++;
                            sumw_(p) += wv;
                            sumww_(p) += wv*wv;
                            actually = true;
                        }
                    }
                    if (actually) {
                        activeAntennas_[icorr].insert(iant);
                        goodRows++;
                    }
                }
            }                
            if (DEVDEBUG && 0) {
                cerr << "flagged " << flagged << endl;
                cerr << "flagSlice " << flagSlice << endl
                     << "fl.shape() " << fl.shape() << endl
                     << "Vpad_.shape() " << Vpad_.shape() << endl
                     << "flagged.shape() " << flagged.shape() << endl
                     << "sl1 " << sl1 << endl;
            }
        }
    }
    if (DEVDEBUG) {
        cerr << "In DelayRateFFTConcat::DelayRateFFTConcat " << endl;
        printActive();
        cerr << "sumw_ " << sumw_ << endl;
        cerr << "Constructed a DelayRateFFTConcat object." << endl;
        cerr << "totalRows " << totalRows << endl;
        cerr << "goodRows " << goodRows << endl;
    }
    
}

DelayRateFFTConcat::DelayRateFFTConcat(Array<Complex>& data, Int nPadFactor, Float f0, Float df, Float dt, SDBList& s,
                           Array<Double>& delayWindow, Array<Double>& rateWindow) :
    DelayRateFFT(s, 0, delayWindow, rateWindow),
    gm_(s),
    nPadFactor_(nPadFactor),
    Vpad_() {
    dt_ = dt;
    f0_ = f0;
    df_ = df;
    IPosition shape = data.shape();
    nCorr_ = shape(0);
    nElem_ = shape(1);
    nt_ = shape(2);
    nChan_ = shape(3);
    nPadT_ = nPadFactor_*nt_;
    nPadChan_ = nPadFactor_*nChan_;
    IPosition paddedDataSize(4, nCorr_, nElem_, nPadT_, nPadChan_);
    Vpad_.resize(paddedDataSize);
    
    IPosition start(4, 0, 0, 0, 0);
    IPosition stop(4, nCorr_,  nElem_, nt_, nChan_);
    IPosition stride(4, 1, 1, 1, 1);
    Slicer sl1(start, stop, stride, Slicer::endIsLength);
    Vpad_(sl1) = data;

    unitize(Vpad_);

}

void
DelayRateFFTConcat::FFT() {
    if (DEVDEBUG) {
        cerr << "DelayRateFFTConcat::FFT()" << endl;
    }
    // Axes are 0: correlation (i.e., hand of polarization), 1: antenna, 2: time, 3: channel
    Vector<Bool> ax(4, false);
    ax(2) = true;
    ax(3) = true;
    // Also copied from DelayFFT in KJones.
    // we make a copy to FFT in place.
    if (DEVDEBUG) {
        cerr << "Vpad_.shape() " << Vpad_.shape() << endl;
    }
    ArrayLattice<Complex> c(Vpad_);
    LatticeFFT::cfft0(c, ax, true);
    if (DEVDEBUG) {
        cerr << "FFT transformed" << endl;
    }
}

std::pair<Bool, Float>
DelayRateFFTConcat::xinterp(Float alo, Float amax, Float ahi) {
    if (amax > 0.0 && alo == amax && amax == ahi)
        return std::make_pair(true, 0.5);

    Float denom(alo-2.0*amax+ahi);
    Bool cond = amax>0.0 && abs(denom)>0.0 ;
    Float fpk = cond ? 0.5-(ahi-amax)/denom : 0.0;
    return std::make_pair(cond, fpk);
}
    
void
DelayRateFFTConcat::searchPeak() {
    if (DEVDEBUG) {
        cerr << "DelayRateFFTConcat::searchPeak()" << endl;
    }

    // Recall param_ -> [phase, delay, rate] for each correlation
    param_.resize(3*nCorr_, nElem_); // This might be better done elsewhere.
    param_.set(0.0);
    flag_.resize(3*nCorr_, nElem_);
    flag_.set(true);  // all flagged initially
    if (DEVDEBUG) {
        cerr << "nt_ " << nt_ << " nPadChan_ " << nPadChan_ << endl;
        cerr << "Vpad_.shape() " << Vpad_.shape() << endl;
        cerr << "delayWindow_ " << delayWindow_ << endl;

    }
    
    for (Int icorr=0; icorr<nCorr_; ++icorr) {
        flag_(icorr*3 + 0, refant()) = false; 
        flag_(icorr*3 + 1, refant()) = false;
        flag_(icorr*3 + 2, refant()) = false;
        for (Int ielem=0; ielem<nElem_; ++ielem) {
            if (ielem==refant()) {
                continue;
            }
            // NB: Time, Channel
            // And once again we fail at slicing
            IPosition start(4, icorr, ielem,      0,         0);
            IPosition stop(4,     1,     1, nPadT_, nPadChan_);
            IPosition step(4,     1,     1,       1,        1);
            Slicer sl(start, stop, step, Slicer::endIsLength);
            Matrix<Complex> aS = Vpad_(sl).nonDegenerate();
            Int sgn = (ielem < refant()) ? 1 : -1;

            // Below is the gory details for turning delay window into index range
            Double bw = Float(nPadChan_)*df_;
            Double d0 = sgn*delayWindow_(IPosition(1, 0));
            Double d1 = sgn*delayWindow_(IPosition(1, 1));
            if (d0 > d1) std::swap(d0, d1);
            d0 = max(d0, -0.5/df_);
            d1 = min(d1, (0.5-1/nPadChan_)/df_);

            // It's simpler to keep the ranges as signed integers and
            // handle the wrapping of the FFT in the loop over
            // indices. Recall that the FFT result returned has indices
            // that run from 0 to nPadChan_/2 -1 and then from
            // -nPadChan/2 to -1, so far as our delay is concerned.
            Int i0 = bw*d0;
            Int i1 = bw*d1;
            if (i1==i0) i1++;
            // Now for the gory details for turning rate window into index range
            Double width = nPadT_*dt_*1e9*f0_;
            Double r0 = sgn*rateWindow_(IPosition(1,0));
            Double r1 = sgn*rateWindow_(IPosition(1,1));
            if (r0 > r1) std::swap(r0, r1);
            r0 = max(r0, -0.5/(dt_*1e9*f0_));
            r1 = min(r1, (0.5 - 1/nPadT_)/(dt_*1e9*f0_));
            
            Int j0 = width*r0;
            Int j1 = width*r1;
            if (j1==j0) j1++;
            if (DEVDEBUG) {
                cerr << "Checking the windows for delay and rate search." << endl;
                cerr << "bw " << bw << endl;
                cerr << "d0 " << d0 << " d1 " << d1 << endl;
                cerr << "i0 " << i0 << " i1 " << i1 << endl; 
                cerr << "r0 " << r0 << " r1 " << r1 << endl;
                cerr << "j0 " << j0 << " j1 " << j1 << endl; 
            }
            Matrix<Float> amp(amplitude(aS));
            Int ipkch(0);
            Int ipkt(0);
            Float amax(-1.0);
            // Unlike KJones we have to iterate in time too
            for (Int itime0=j0; itime0 != j1; itime0++) {
                Int itime = (itime0 < 0) ? itime0 + nPadT_ : itime0;
                for (Int ich0=i0; ich0 != i1; ich0++) {
                    Int ich = (ich0 < 0) ? ich0 + nPadChan_ : ich0;
                    // cerr << "Gridpoint " << itime << ", " << ich << "->" << amp(itime, ich) << endl;
                    if (amp(itime, ich) > amax) {
                        ipkch = ich;
                        ipkt  = itime;
                        amax=amp(itime, ich);
                    }
                }
            }
            // Finished grovelling. Now we have the location of the
            // maximum amplitude.
            Float alo_ch = amp(ipkt, (ipkch > 0) ? ipkch-1 : nPadChan_-1);
            Float ahi_ch = amp(ipkt, ipkch<(nPadChan_-1) ? ipkch+1 : 0);
            std::pair<Bool, Float> maybeFpkch = xinterp(alo_ch, amax, ahi_ch);
            // We handle wrapping while looking for neighbours
            Float alo_t = amp(ipkt > 0 ? ipkt-1 : nPadT_ -1,     ipkch);
            Float ahi_t = amp(ipkt < (nPadT_ -1) ? ipkt+1 : 0,   ipkch);
            if (DEVDEBUG) {
                cerr << "For element " << ielem << endl;
                cerr << "In channel dimension ipkch " << ipkch << " alo " << alo_ch
                     << " amax " << amax << " ahi " << ahi_ch << endl;
                cerr << "In time dimension ipkt " << ipkt << " alo " << alo_t
                     << " amax " << amax << " ahi " << ahi_t << endl;
            }
            std::pair<Bool, Float> maybeFpkt = xinterp(alo_t, amax, ahi_t);

            if (maybeFpkch.first and maybeFpkt.first) {
                // Phase
                Complex c = aS(ipkt, ipkch);
                Float phase = arg(c);
                param_(icorr*3 + 0, ielem) = sgn*phase;
                Float delay = (ipkch)/Float(nPadChan_);
                if (delay > 0.5) delay -= 1.0;           // fold
                delay /= df_;                           // nsec
                param_(icorr*3 + 1, ielem) = sgn*delay; //
                Double rate = (ipkt)/Float(nPadT_);
                if (rate > 0.5) rate -= 1.0;
                Double rate0 = rate/dt_;
                Double rate1 = rate0/(1e9 * f0_); 

                param_(icorr*3 + 2, ielem) = Float(sgn*rate1); 
                if (DEVDEBUG) {
                    cerr << "maybeFpkch.second=" << maybeFpkch.second
                         << ", df_= " << df_ 
                         << " fpkch=" << (ipkch + maybeFpkch.second) << endl;
                    cerr << " maybeFpkt.second=" << maybeFpkt.second
                         << " rate0=" << rate
                         << " 1e9 * f0_=" << 1e9 * f0_ 
                         << ", dt_=" << dt_
                         << " fpkt=" << (ipkt + maybeFpkt.second) << endl;
                        
                }
                if (DEVDEBUG) {
                    cerr << "Found peak for element " << ielem << " correlation " << icorr
                         << " ipkt=" << ipkt << "/" << nPadT_ << ", ipkch=" << ipkch << "/" << nPadChan_
                         << " peak=" << amax 
                         << "; delay " << delay << ", rate " << rate
                         << ", phase " << arg(c) << " sign= " << sgn << endl;
                }
                // Set 3 flags.
                flag_(icorr*3 + 0, ielem)=false; 
                flag_(icorr*3 + 1, ielem)=false;
                flag_(icorr*3 + 2, ielem)=false;
            }
            else {
                if (DEVDEBUG) {
                    cerr << "No peak in 2D FFT for element " << ielem << " correlation " << icorr << endl;
                }
            }
        }
    }
}


Float
DelayRateFFTConcat::snr(Int icorr, Int ielem, Float delay, Float rate) {
    if (DEVDEBUG) {
        cerr << "DelayRateFFTConcat::snr"<< endl;
    }
    // We calculate a signal-to-noise ration for the 2D FFT fringefit
    // using a formula transcribed from AIPS FRING.
    //
    // Have to convert delay and rate back into indices on the padded 2D grid.
    Int sgn = (ielem < refant()) ? 1 : -1;
    delay *= sgn*df_;
    if (delay < 0.0) delay += 1;
    Int ichan = Int(delay*nPadChan_ + 0.5); 
    if (ichan == nPadChan_) ichan = 0;
        
    rate *= sgn*1e9 * f0_;
    rate *= dt_;
    if (rate < 0.0) rate += 1;
    Int itime = Int(rate*nPadT_ + 0.5);
    if (itime == nPadT_) itime = 0;
    // What about flags? If the datapoint closest to the computed
    // delay and rate values is flagged we probably shouldn't use
    // it, but what *should* we use?
    IPosition ipos(4, icorr, ielem, itime, ichan);
    IPosition p(2, icorr, ielem);
    Complex v = Vpad_(ipos);
    Float peak = abs(v);
    if (peak > 0.999*sumw_(p)) peak=0.999*sumw_(p);
    // xcount is number of data points for baseline to ielem
    // sumw is sum of weights,
    // sumww is sum of squares of weights

    Float cwt;
    if (fabs(sumw_(p))<FLT_EPSILON) {
        cwt = 0;
        if (DEVDEBUG) {
            cerr << "Correlation " << icorr << " antenna " << ielem << ": sum of weights is zero." << endl;
        }
    } else {
        Float x = C::pi/2*peak/sumw_(p);
        // The magic numbers in the following formula are from AIPS FRING
        cwt = (pow(tan(x), 1.163) * sqrt(sumw_(p)/sqrt(sumww_(p)/xcount_(p))));
        if (DEVDEBUG) {
            cerr << "Correlation " << icorr << " antenna " << ielem << " ipos " << ipos
                 << " peak=" << peak << "; xang=" << x << "; xcount=" << xcount_(p) << "; sumw=" << sumw_(p) << "; sumww=" << sumww_(p)
                 << " snr " << cwt << endl;
        }
    }
    return cwt;
}
    

} // end namespace
