/**
   \file lineshapes.cpp
   Bojan Nikolic <bn204@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   
   Initial version February 2008
   Maintained by ESO since 2013.
   Renamed lineshapes.cc 2023

*/

#include "lineshapes.h"
#include "lineparams.h"

namespace LibAIR2 {

  double CGrossLine(double f,
		    const CLineParams & cp)
  {
    return GrossLine( f,
		      cp.f0, cp.gamma, cp.S);
  }

}


