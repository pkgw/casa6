/**
   \file models_basic.cpp
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   
   Initial version February 2008
   Revised 2009
   Maintained by ESO since 2013.

   Renamed models_basic.cc 2023

*/

#include <stdexcept>

#include "models_basic.h"

#include "radiometermeasure.h"
#include "slice.h"
#include "columns.h"
#include "lineparams.h"
#include "basicphys.h"
#include "partitionsum.h"
#include "partitionsum_testdata.h"
#include "rtranfer.h"
#include "layers.h"

namespace LibAIR2 {

  std::shared_ptr<Radiometer> SwitchRadiometer(RadiometerT r)
  {
    Radiometer *res;
    switch (r)
    {
    case ALMAProd:
      res=MkFullALMAWVR();
      break;
    case ALMADickeProto:
      res=MkFullDickeProtoWVR();
      break;
    case IRAM22GHz:
      res=MkIRAM22();
      break;
    default:
      throw std::runtime_error("Unknown radiometer type");
    }
    return std::shared_ptr<Radiometer>(res);
  }
}


