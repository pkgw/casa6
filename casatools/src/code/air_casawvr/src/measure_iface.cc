/**
   \file measure_iface.cpp
   Renamed measure_iface.cc 2023
   
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>

*/

#include "measure_iface.h"
#include "model_iface.h"
#include "numalgo.h"

namespace LibAIR2 {

  ALMAMeasure::ALMAMeasure(WVRAtmoQuantModel &model):
    model(model)
  {
  }
  
  void ALMAMeasure::AddParams ( std::vector< Minim::DParamCtr > &pars )
  {
    model.AddParams(pars);
  }

  NormalNoise::NormalNoise(size_t n):
    thermNoise(n)
  {
  }
  

  AbsNormMeasure::AbsNormMeasure(WVRAtmoQuantModel &model):
    ALMAMeasure(model),
    NormalNoise(4),
    obs(4)
  {
  }

  void AbsNormMeasure::modelObs(void)
  {
    model.eval(obs);
  }

  double AbsNormMeasure::lLikely (void) const
  {
    std::vector<double> res;
    model.eval(res);
    
    return GaussError( obs.begin(), res, thermNoise);
    
  }

}

