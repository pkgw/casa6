/**
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   Initial version November 2009
   Maintained by ESO since 2013.

   \file almaabs_i.hpp
   Renamed almaabs_i.h 2023


*/
#ifndef _LIBAIR_ALMAABS_I_HPP__
#define _LIBAIR_ALMAABS_I_HPP__

#include <vector>
#include <list>

#include <memory>

#include "bnmin1/src/nestedsampler.h"
#include "bnmin1/src/priors.h"

#include "almaabs.h"
#include "../model_iface.h"
#include "../dipmodel_iface.h"
#include "../measure_iface.h"

#include "almaresults.h"


namespace LibAIR2 {

  /// Structures to represent likelihood of a measurement for an
  /// absolute retrieval from ALMA data
  struct iALMAAbsRetLL
  {
    /// Assume thermal noise 
    static const double thermNoise;

    /// Model of the atmosphere before taking into account elevation
    /// of observation
    CouplingModel *cm;

    /// Model of the atmosphere after taking into account elevation
    PPDipModel m;

    /// Representation of the measured values and errors
    AbsNormMeasure *ll;

    iALMAAbsRetLL(const std::vector<double> &TObs,
		  double el,
		  const ALMAWVRCharacter &WVRChar);

  };

  class iALMAAbsRet
  {

  public:

    iALMAAbsRetLL ls;

    /// The posterior 
    std::list<Minim::WPPoint> post;

    /// Representation of the likelihood and priors 
    Minim::IndependentFlatPriors pll;

    /// Evidence value
    double evidence;

    /// The nested sampler
    std::unique_ptr<Minim::NestedS> ns;
    

    /// Number of points in the live set
    static const size_t n_ss;

    iALMAAbsRet(const std::vector<double> &TObs,
		double el,
		const ALMAWVRCharacter &WVRChar);

    bool sample(void); // returns false if evidence is zero

    // -------------- Retrieval of results ------------------
    
    /** Get the important model parameters and estimated errors
     */
    void  g_Pars(ALMAResBase &r);

    /** Get the inferred phase-correction coefficients and estimated
	errors
    */
    void  g_Coeffs(ALMAResBase &r);
    
  };


}



#endif
