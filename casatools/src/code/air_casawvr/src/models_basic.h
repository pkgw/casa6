/**
   Bojan Nikolic <bn204@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   Initial version February 2008
   Revised 2009
   Maintained by ESO since 2013.

   This file is part of LibAIR and is licensed under GNU Public
   License Version 2

   \file models_basic.hpp
   Renamed models_basic.h 2023

   Models which expose conveniently a limited set of variable
   parameters and encapsulate other complexities.
*/

#ifndef _LIBAIR_MODELS_BASIC_HPP__
#define _LIBAIR_MODELS_BASIC_HPP__

#include <vector>
#include <memory>

#include "bnmin1/src/minimmodel.h"
#include "singlelayerwater.h"

#include "model_iface.h"


namespace LibAIR2 {

  // Forward decleration
  class Radiometer;


  /** \brief Create radiometer object given model design

      \returns the correct radiometer according to type supplied
   */
  std::shared_ptr<Radiometer> SwitchRadiometer(RadiometerT r);


}

#endif
