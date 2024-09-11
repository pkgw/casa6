/** 
    \file twoerrline_imp.hxx

    Bojan Nikolic <bojan@bnikolic.co.uk>, <b.nikolic@mrao.cam.ac.uk>
    Initial version 2009    
    
    Maximum likelihood tools for fitting a line to two noisy
    coordinates

    Renamed to twoerrline_ml.h 2023

*/
#ifndef _BNMIN1_TWOERRLINE_ML_HXX__
#define _BNMIN1_TWOERRLINE_ML_HXX__

#include <casacore/casa/Arrays/MatrixMath.h>

#include "linemodel.h"
#include "minimmodel.h"

namespace Minim {

  /**
     \brief Fit straight-line model to data with noise in both coordinates
     
   */
  class LineTwoErrML:
    public LineMod
  {

    casacore::Vector<double> xobs;
    casacore::Vector<double> yobs;

  public:

    // ---------- Public data  -----------------------    

    const double sigmax;

    const double sigmay;

    const size_t nobs;

    // ------------ Construction / Destruction ----------------
    
    LineTwoErrML(const std::vector<double> &xvals,
		 const std::vector<double> &yvals,
		 double sigmax,
		 double sigmay);

    // ------------- Public interface -------------------------

    /// Evalute the residuals from the current line
    void residuals(casacore::Vector<double> &res) const;

    // -------------------- Inherited functions ---------------
    double lLikely(void) const;
    void lGrd(std::vector< double > &res) const;
  };


  /**
     \brief Fit straight-line model to data with noise in both
     coordinates, using the Levengberg-Marquand algorithm
     
   */
  class LineTwoErr_LavMarq:
    public Minimisable
  {
    LineTwoErrML m;
  public:
    
    LineTwoErr_LavMarq(const std::vector<double> &xvals,
		       const std::vector<double> &yvals,
		       double sigmax,
		       double sigmay);

    // -------------------- Inherited functions ---------------
    void AddParams(std::vector< Minim::DParamCtr > &pars);    
    virtual void  residuals (std::vector<double> &res) const;    
    virtual unsigned nres (void) const;     

  };


}

#endif
