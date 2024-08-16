// Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
// Initial version July 2009
// Maintained by ESO since 2013.

/** \file model_make.cpp

    Renamed model_make.cc 2023

 */

#include "model_make.h"
#include "model_water.h"
#include "singlelayerwater.h"
#include "models_basic.h"
#include "slice.h"
#include "layers.h"
#include "rtranfer.h"
#include "radiometermeasure.h"

namespace LibAIR2 {

  LibAIR2::WaterModel<ISingleLayerWater> *
  mkSingleLayerWater(RadiometerT radiot, 
		     PartitionTreatment t,
		     Continuum c,
		     double PDrop)
  {
    std::shared_ptr<Radiometer> r(SwitchRadiometer(radiot));

    std::shared_ptr<ISingleLayerWater> sl(new ISingleLayerWater (r->getFGrid(),
								 WaterData::L183,
								 t,
								 c,
								 PDrop));
    return new LibAIR2::WaterModel<ISingleLayerWater> (r,
						       sl);
    
  }

  LibAIR2::WaterModel<ISingleLayerWater> *
  mkSingleLayerWater(const ALMAWVRCharacter &ac, 
		     PartitionTreatment t,
		     Continuum c,
		     double PDrop)
  {
    std::shared_ptr<Radiometer> r(MkALMAWVR(ac));

    std::shared_ptr<ISingleLayerWater> sl(new ISingleLayerWater (r->getFGrid(),
								 WaterData::L183,
								 t,
								 c,
								 PDrop));
    return new LibAIR2::WaterModel<ISingleLayerWater> (r,
						       sl);

  }

  LibAIR2::WaterModel<ISingleLayerWater> *
  mkSimpleOffset(double cf,
		 double bw)
  {
    std::shared_ptr<Radiometer> r(MkALMAWVR_offset(cf,bw));

    std::shared_ptr<ISingleLayerWater> sl(new ISingleLayerWater (r->getFGrid(),
								 WaterData::L183,
								 PartTable,
								 AirCont,
								 0));
    return new LibAIR2::WaterModel<ISingleLayerWater> (r,
						       sl);

  }

  LibAIR2::WaterModel<ICloudyWater> *
  mkCloudy(RadiometerT radiot, 
	   PartitionTreatment t,
	   Continuum c,
	   double PDrop)
  {
    std::shared_ptr<Radiometer> r(SwitchRadiometer(radiot));

    std::shared_ptr<ICloudyWater> sl(new ICloudyWater (r->getFGrid(),
						       WaterData::L183,
						       t,
						       c,
						       PDrop));
    return new LibAIR2::WaterModel<ICloudyWater> (r,
						  sl);

  }

  LibAIR2::WaterModel<ICloudyWater> *
  mkCloudy(const ALMAWVRCharacter &ac, 
	   PartitionTreatment t,
	   Continuum c,
	   double PDrop)
  {
    std::shared_ptr<Radiometer> r(MkALMAWVR(ac));

    std::shared_ptr<ICloudyWater> sl(new ICloudyWater (r->getFGrid(),
						       WaterData::L183,
						       t,
						       c,
						       PDrop));
    return new LibAIR2::WaterModel<ICloudyWater> (r,
						  sl);

  }
  

}



