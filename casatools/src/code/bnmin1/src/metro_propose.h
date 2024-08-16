/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2008

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2

   \file metro_propose.hxx
   Renamed to metro_propose.h

   Calculation of proposal points for Metropolis sampling
*/
#ifndef _BNMIN1_METRO_PROPOSAL_HXX__
#define _BNMIN1_METRO_PROPOSAL_HXX__

#include <vector>

#include <random>

namespace Minim {

  /** \brief A helper class which handles the random number parts of
      the Metropolis algorithm

      The two objectives of this class functions are to generate the
      proposal points [via the function displace()] and to generate
      uniform random numbers for accepting less likely points

      It encapsulates a random number generator 
   */
  class MetroPropose {

  protected:
    
    /// Standard distributions of the proposal points
    std::vector<double> sigmas;

    // Stuff for random numbers
    typedef std::mt19937  base_generator_type;
    base_generator_type generator;
    std::normal_distribution<double> norm_dist;
    std::uniform_real_distribution<double> uni_dist;
        
    
  public:
    
    // ---------- Construction / Destruction --------------

    /**
       \param sigmas The standard distribution to be used for each of
       the parameters that is being varied in the model
     */
    MetroPropose(const std::vector<double> & sigmas,
		 unsigned seed=0);
    
    virtual ~MetroPropose();

    // ---------- Public interface --------------------------

    /// Normal distribution for creating proposal points
    double norm();

    /// Uniform distribution for accepting according to probability
    double uni();
    
    /**
       Generate a proposad point by displacing the point x
    */
    virtual void displace(std::vector<double> &x);

    /**
       \returns number of parameters to propose for
    */
    size_t nPars(void);

    /**
      Generate a uniform number between zero and one for acceptance
      probability
    */
    double raccept(void)
    {
      return uni();
    }

    /** \brief Scale the variance of the of the proposal distributions
     */
    void scaleSigma(double c);
  
  };

  /** \brief Proposes each parameter in turn

     Intead of displacing the entire vector at once, each parameter is
     displaced in turn. See MetropolisMCMC::Options::Sequence.
   */
  class MetroProposeSeq :
    public MetroPropose
  {

    size_t count;

    const size_t n;

  public:

    // ---------- Construction / Destruction --------------

    MetroProposeSeq(const std::vector<double> & sigmas,
		    unsigned seed=0);

    // ---------- Public interface --------------------------
    virtual void displace( std::vector<double> &x);

  };

}

#endif
