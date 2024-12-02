/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2008

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2

   \file metro_propose.cxx
   Renamed to metro_propose.cc

*/

#include "metro_propose.h"

namespace Minim {

  MetroPropose::MetroPropose(const std::vector<double> & sigmas,
			     unsigned seed):
    sigmas(sigmas),
    norm_dist(0.,1.),
    uni_dist(0.,1.)
  {
    if(seed != 0){
      generator.seed(seed);
    }
  }

  MetroPropose::~MetroPropose()
  {
  }

  double MetroPropose::norm(){
    return norm_dist(generator);
  }

  double MetroPropose::uni(){
    return uni_dist(generator);
  }
  
  void MetroPropose::displace( std::vector<double> &x)
  {
    for (size_t i =0 ; i < sigmas.size() ; ++i)
      x[i] += sigmas[i]* norm();
  }

  size_t MetroPropose::nPars(void)
  {
    return sigmas.size();
  }

  void MetroPropose::scaleSigma(double c)
  {
    for(size_t i=0; i<sigmas.size(); ++i)
      sigmas[i] *= c;
  }

  MetroProposeSeq::MetroProposeSeq(const std::vector<double> & sigmas,
				   unsigned seed):
    MetroPropose(sigmas, seed),
    count(0),
    n(sigmas.size())
  {
  }

  void MetroProposeSeq::displace( std::vector<double> &x)
  {
    const size_t i = count  % n;
    x[i] += sigmas[i]* norm();
    ++count;
  }
    

}

