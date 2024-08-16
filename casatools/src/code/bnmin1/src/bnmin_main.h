/**
   Bojan Nikolic <bojan@bnikolic.co.uk> 
   Initial version 2008

   This file is part of BNMin1 and is licensed under GNU General
   Public License version 2.

   \file bnmin_main.hxx
   Renamed to bnmin_main.h 2023.

   The main include file for BNMin1 Library
   
   \mainpage A simple minimisation / inference library

*/

#include <string>
#include <stdexcept>

#ifndef __BNMIN_BNMIN_MAIN_HPP__
#define __BNMIN_BNMIN_MAIN_HPP__

/* Define to the full name of this package. */
#define PACKAGE_NAME "BNMIN1"
/* Define to the full name and version of this package. */
#define PACKAGE_STRING "BNMIN1 1.11"
/* Define to the version of this package. */
#define PACKAGE_VERSION "1.11"


namespace Minim {


  const char * version(void);

  /** \brief Base class for run-time errors within the library 
   */
  class BaseErr:
    public std::runtime_error
  {
  public:
    BaseErr(const std::string &s);

  };

}
#endif

