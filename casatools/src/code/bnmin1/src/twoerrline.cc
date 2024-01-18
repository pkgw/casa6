/**

   Bojan Nikolic <bojan@bnikolic.co.uk>, <b.nikolic@mrao.cam.ac.uk>
   Initial version 2009

   \file twoerrline.cxx
   Renamed to twoerrline.cc 2023

   Fit lines to data which have errors in both coordinates
*/

#include <memory>
#include "twoerrline.h"
#include "twoerrline_ml.h"
#include "gradientminim.h"
#include "lmmin.h"
#include "lmminutils.h"

namespace Minim {

  void LFit_MaxL(const std::vector<double> &xvals,
		 const std::vector<double> &yvals,
		 double sigma_x,
		 double sigma_y,
		 LineFit &res)
  {
    // We dont get the covariance matrix using the BFGS2 algorithm
    res.cv[0]=res.cv[3]=-1;

    Minim::LineTwoErrML lml(xvals, 
			    yvals,
			    sigma_x, 
			    sigma_y);
    lml.a=res.a;
    lml.b=res.b;

    Minim::BFGS2Minim m(lml);
    m.solve();
    
    res.a=lml.a;
    res.b=lml.b;
    
  }

  void LFit_LM(const std::vector<double> &xvals,
	       const std::vector<double> &yvals,
	       double sigma_x,
	       double sigma_y,
	       LineFit &res)
  {
    Minim::LineTwoErr_LavMarq lml(xvals, 
				  yvals,
				  sigma_x,
				  sigma_y);
    Minim::LMMin m(lml);

    m.getbyname("a")->setp(res.a);
    m.getbyname("b")->setp(res.b);

    m.solve();
    
    res.a=m.getbyname("a")->getp();
    res.b=m.getbyname("b")->getp();
    
    std::unique_ptr<std::vector<double> > cv(CVMatrix(m,1e-5));

    std::copy(cv->begin(),
	      cv->end(),
	      res.cv.begin());
    

  }

}
