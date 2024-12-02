/**
   
   Bojan Nikolic <bojan@bnikolic.co.uk>, <b.nikolic@mrao.cam.ac.uk>

   \file paramalgo.cxx
   Renamed to paramalgo.cc 2023
   
*/

#include "paramalgo.h"

namespace Minim {

  DParamCtr * findName(std::vector<DParamCtr> & parv,
		       const std::string & pname)
  {
    for ( std::vector<DParamCtr>::iterator i ( parv.begin() ) ;
	  i < parv.end() ;
	  ++i ) 
      if (i->name == pname ) return &(*i);
    
    return NULL;
  }

}

