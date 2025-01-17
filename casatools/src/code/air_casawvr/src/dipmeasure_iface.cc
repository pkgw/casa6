/**
   \file dipmeasure_iface.cpp

   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk> 

   Renamed dipmeasure_iface.cc 2023
  
*/

#include "dipmeasure_iface.h"
#include "dipmodel_iface.h"
#include "numalgo.h"

namespace LibAIR2{
  
  DipNormMeasure::DipNormMeasure(PPDipModel & model):
    ALMAMeasure(model),
    NormalNoise(4),
    _model(model)
  {
  }

  void DipNormMeasure::addObs(double za,
			      const std::vector<double> & skyTb)
  {
    obs.push_back( std::make_tuple(za, skyTb));
  }

  double DipNormMeasure::lLikely (void) const
  {
    std::vector<double> scratch;
    double res=0;
    for (size_t i = 0 ; i < obs.size() ; ++i)
    {
      double za;
      obs_t  skyT;
      std::tie(za, skyT) = obs[i];

      _model.setZA(za);
      _model.eval(scratch);

      res+= GaussError( skyT.begin(), scratch, thermNoise);
    }

    return res;
  }
  
  

}

