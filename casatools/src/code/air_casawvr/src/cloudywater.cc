// Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
// Initial version June 2009
// Maintained by ESO since 2013.

/** \file cloudywater.cpp
    Renamed cloudywater.cc 2023
    
 */

#include "cloudywater.h"

#include "layers.h"
#include "slice.h"
#include "rtranfer.h"

namespace LibAIR2 {

  ICloudyWater::ICloudyWater(const std::vector<double> &fgrid,
			     LibAIR2::WaterData::Lines l,
			     PartitionTreatment t,
			     Continuum c,
			     double PDrop):
    ISingleLayerWater(fgrid,
		      l,
		      t,
		      c,
		      PDrop),
    tau183(0),
    CT(270),
    contslice(CT,0),
    sr(contslice,
       fgrid),
    col(0, 183.3)
  {
    contslice.AddColumn(col);
  }

  const std::vector<double> &   ICloudyWater::TbGrid(void)
  {
    updatePars();
    sr.UpdateI(ISingleLayerWater::getBckg());
    return ISingleLayerWater::TbGrid(sr);
  }

  void ICloudyWater::updatePars(void)
  {
    col.setN(tau183);
    contslice.setT(CT);
  }

}
