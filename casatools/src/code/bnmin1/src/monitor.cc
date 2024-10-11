/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2003

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2

   \file monitor.cxx
   Renamed to monitor.cc 2023

*/

#include "monitor.h"

#include "minim.h"
#include "minimio.h"

#include <iostream>

namespace Minim {


  void ChiSqMonitor::Iter ( Minimiser * m)
  {
    std::cerr<<"Chi-squared: " << m->ChiSquared() <<std::endl;
  }

  

  void  ParsMonitor::Iter ( Minimiser * m)
  {
    PrettyPrint ( *m );
  }

}



