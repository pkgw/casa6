/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2008

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2.

   \file bnmin_main.cxx
   Renamed to bnmin_main.cc 2023.

*/


#include "bnmin_main.h"

namespace Minim {

  const char * version(void)
  {
    return PACKAGE_VERSION;
  }

  BaseErr::BaseErr(const std::string &s):
    std::runtime_error(s)
  {
  }


}


