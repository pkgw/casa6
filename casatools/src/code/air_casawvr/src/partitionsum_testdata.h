/**
   Bojan Nikolic <bn204@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   Initial version February 2008
   Maintained by ESO since 2013.

   This file is part of LibAIR and is licensed under GNU Public
   License Version 2

   \file partitionsum_testdata.hpp
   Renamed partitionsum_testdata.h 2023

   Some data for testing/basic usage of partition tables
*/

#ifndef __LIBAIR_PARTITIONSUM_TESTDATA_HPP__
#define __LIBAIR_PARTITIONSUM_TESTDATA_HPP__

#include "partitionsum.h"

namespace LibAIR2 {

  /**
     Return the raw partition table for h2o
   */
  casacore::Matrix<double>* getH2ORawTable(void);


}

#endif
