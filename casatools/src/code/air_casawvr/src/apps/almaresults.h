/**
   Bojan Nikolic <b.nikolic@mrao.cam.ac.uk>, <bojan@bnikolic.co.uk>
   Initial version November 2009
   Maintained by ESO since 2013.

   \file almaresults.hpp

*/
#ifndef _LIBAIR_APPS_ALMARESULT_HPP__
#define _LIBAIR_APPS_ALMARESULT_HPP__

#include <vector>
#include <iosfwd>
#include <list>

#include "alma_datastruct.h"

namespace LibAIR2 {

  /** \brief Definition of quantities to be typically retrieved from
      ALMA analysis
  */
  struct ALMAResBase: 
    public ALMARes_Basic 
  {

    /**
       
     */
    ALMAResBase(void);

    virtual ~ALMAResBase();
    
    /**\brief Output the results separated by tabs and on a single line 

     */
    void print_str_inline(std::ostream &os) const;
    
    
  };

  struct ALMAResBaseList // simple replacement for what we need from boost::ptr_list
  {

    ALMAResBaseList(void);

    virtual ~ALMAResBaseList();
    
    std::list<ALMAResBase*> ptr_list;
  };
  
  std::ostream &operator<<(std::ostream &os,
			   const ALMAResBase &r);

  std::ostream &operator<<(std::ostream &os,
			   const std::list<ALMAResBase*> &i);

  struct ALMAContRes:
    public ALMAResBase {
    
    /** \brief Estimate of non-water vapour continuum opacity at 183
       ghz
    */
    double tau183;

    double tau183_err;

    std::ostream &header_inline(std::ostream &os) const;
    virtual void print_str_inline(std::ostream &os);
  };

  std::ostream &operator<<(std::ostream &os,
			   const ALMAContRes &r);
    


}

#endif


