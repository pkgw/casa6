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
#include <sstream>

#include <stdio.h> 
#include <stdlib.h>
#include <string>

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
    if (not ifs.good()){
      throw std::runtime_error(std::string("Could not open dispersion table ")+fname);
    }
    std::string scratch;

    while(ifs.good()){
      std::getline(ifs, scratch);
      //std::cerr << scratch << " - interpreted as:"; 
      if (scratch.size() < 5){
	//std::cerr << "(nothing)" << std:: endl;
	continue;
      }

      std::stringstream ss(scratch);
      double first, second;
      std::string sep;
      if(ss >> first){
	if(ss >> sep){
	  if(ss >> second){
	    dt.insert(dt.end(),
		      std::pair<double, double>(first, second));
	    //std::cerr << "(double)  " << first << " ,  " << second << " " << std::endl;
	  }
	  else{
	    std::cerr<<" Reading " << fname << ": could not interpret third part of " << scratch <<std::endl;
	  }
	}
	else{
	  std::cerr<<" Reading " << fname << ": could not interpret separator in " << scratch <<std::endl;
	}
      }
      else{
	std::cerr<<" Reading " << fname << ": could not interpret first part of " << scratch <<std::endl;
      } 
    }

    return;
  }

}
