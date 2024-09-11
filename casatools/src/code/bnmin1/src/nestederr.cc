/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2009

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2

   \file nestederr.cxx
   Renamed to nestederr.cc 2023


*/

#include "nestederr.h"

namespace Minim {

  NestedSmallStart::NestedSmallStart(const std::list<MCPoint> & /*ss*/):
    BaseErr("Number of points in the starting set is less than two")
  {
  }


}


