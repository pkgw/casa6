/**
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   Initial version August 2010.
   Maintained by ESO since 2013.

   This file is part of LibAIR and is licensed under GNU Public
   License Version 2

   \file dispersion.cpp
   Renamed dispersion.cc 2023

*/
#include <iostream>
#include <fstream>

#include <stdio.h> 
#include <stdlib.h>
#include <string>
#include <cstring>
#include <casacore/casa/BasicSL/String.h>

#include "dispersion.h"

namespace LibAIR2 {

  double DispersionTab::operator() (double fnu)
  {
    std::pair<double, double> b= *lower_bound(fnu);
    if (b.first==fnu)
      return b.second;

    std::pair<double, double> l= *(--lower_bound(fnu));
    std::pair<double, double> u= *upper_bound(fnu);
    
    const double f=(fnu-l.first)/(u.first-l.first);
    return l.second+ f*(u.second-l.second);
    
  }

  void loadCSV(const char *fname,
	       DispersionTab &dt)
  {
    std::ifstream ifs(fname);
    if (not ifs.good())
    {
      throw std::runtime_error(std::string("Could not open dispersion table ")+fname);
    }
    std::string scratch;

    while(ifs.good())
    {
      std::getline(ifs, scratch);
      if (scratch.size() < 5)
	continue;

      char * pch;
      char *s = new char[scratch.size()+1];
      strcpy( s, scratch.c_str() );
      pch = strtok(s, ",;\"");
      casacore::String first(*pch);
      pch = strtok (NULL, ",;\"");
      casacore::String second(*pch);
      
      casacore::trim(first);
      casacore::trim(second);
      try {
	dt.insert(dt.end(),
		  std::pair<double, double>(casacore::String::toDouble(first),
					    casacore::String::toDouble(second)
					    ));
      }
      catch (const std::bad_cast &bc)
      {
	std::cerr<<"Could not interpret " << first << " and " << second
		 <<std::endl;
      }
      delete [] s;
    }
  }


}



