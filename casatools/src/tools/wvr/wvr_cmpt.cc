#include <iostream>
#include <sstream>
#include <fstream>
#include <numeric>
#include <memory>
#include <tuple>
#include <vector>
#include <string>

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Arrays/ArrayUtil.h>
#include <casacore/ms/MeasurementSets/MeasurementSet.h>
#include <casacore/casa/IO/ArrayIO.h>

#include <air_casawvr/casawvr/mswvrdata.h>
#include <air_casawvr/casawvr/msgaintable.h>
#include <air_casawvr/casawvr/msutils.h>
#include <air_casawvr/casawvr/msspec.h>
#include <air_casawvr/casawvr/msantdata.h>
#include <air_casawvr/src/apps/arraydata.h>
#include <air_casawvr/src/apps/arraygains.h>
#include <air_casawvr/src/apps/almaabs.h>
#include <air_casawvr/src/apps/dtdlcoeffs.h>
#include <air_casawvr/src/apps/almaresults.h>
#include <air_casawvr/src/apps/segmentation.h>
#include <air_casawvr/src/libair_main.h>

#include <asdmstman/AsdmStMan.h>

#include <casacore/casa/Logging/StreamLogSink.h>
#include <casacore/casa/Logging/LogSink.h>
#include <casacore/casa/Logging/LogIO.h>

#include <wvr_cmpt.h>

#include "wvrgcalerrors.h"
#include "wvrgcalfeedback.h"

using casacore::LogIO;

// Utilities ////////////////////////////////////////////////////////////

/** \brief Take a parameter than can be specified as a sequence of
    antenna numbers or names and always return as a sequence of
    antenna numbers.
 */						   
std::vector<size_t> getAntParsV(const casacore::String &s,
				const casacore::MeasurementSet &ms)
{
  using namespace LibAIR2;
  aname_t anames=getAName(ms);
  casacore::Vector<casacore::String> pars = stringToVector(s, ',');
  std::vector<size_t> res;

  for (size_t i=0; i<pars.size(); ++i)
  {
    size_t match;
    std::string thisant = pars[i];
    if (std::accumulate( anames.begin( ),
			 anames.end( ), false,
			 [&](bool acc, const aname_t::value_type &p){
			   bool cmp = p.second == thisant;
			   if ( cmp ) match = p.first;
			   return acc || cmp;
			 }
			 ))
      {
	res.push_back(match);
      }
    else
      {
	// should be an antenna number
	try {
	  int n=casacore::String::toInt(thisant, true);
	  if ( std::accumulate( anames.begin(),
				anames.end(), false,
				[=](bool acc, const aname_t::value_type &p) {
				  return acc || ((int)p.first == n);
				}
				) == false )
	    {
	      throw AntIDError(n,
			       anames);
	    }
	  res.push_back(n);
	}
	catch (casacore::AipsError x)
	  {
	    throw AntIDError(thisant,
			     anames); 
	  }	
      }
  }
  return res;
}

LibAIR2::AntSet getAntPars(const casacore::String &s,
			   const casacore::MeasurementSet &ms)
{
  std::vector<size_t> resv = getAntParsV(s, ms);
  LibAIR2::AntSet res;
  

  for (size_t i=0; i<resv.size(); ++i)
  {
    res.insert(resv[i]);
  }

  return res;
}



/* Determine the nearest n antennas not accepting ones which
   are flagged or have a distance > maxdist_m
*/
 
LibAIR2::AntSetWeight limitedNearestAnt(const LibAIR2::antpos_t &pos,
				       size_t i,
				       const LibAIR2::AntSet &flag,
				       size_t n,
				       double maxdist_m)
{
  LibAIR2::AntSetD dist=LibAIR2::antsDist(pos, i, flag);
  LibAIR2::AntSetWeight res;
    
  double total=0;
  size_t limitedn=0;
  LibAIR2::AntSetD::const_iterator s=dist.begin();
  for (size_t j=0; j<n; ++j)
  {
    if(s!=dist.end() and s->first <= maxdist_m)
    {
      total+=s->first;
      ++s;
      ++limitedn;
    }
  }

  s=dist.begin();
  for (size_t j=0; j<limitedn; ++j)
  {
    res.insert(std::make_pair(s->first/total, s->second));
    ++s;
  }

  return res;

}


/** \brief Flag and interpolate WVR data
 */
void flagInterp(const casacore::MeasurementSet &ms,
		const LibAIR2::AntSet &wvrflag,
		LibAIR2::InterpArrayData &d,
		const double maxdist_m,
		const int minnumants,
		LibAIR2::AntSet &interpImpossibleAnts,
		casacore::LogIO *itsLog)
{

  LibAIR2::antpos_t apos;
  LibAIR2::getAntPos(ms, apos);
  LibAIR2::AntSet wvrflag_s(wvrflag.begin(), 
			   wvrflag.end());

  for(LibAIR2::AntSet::const_iterator i=wvrflag.begin();
      i!=wvrflag.end(); 
      ++i)
  {

    LibAIR2::AntSetWeight near=limitedNearestAnt(apos, 
						 *i, 
						 wvrflag_s, 
						 3,
						 maxdist_m);
    if(near.size()>= static_cast<unsigned int>(minnumants)){
      //LibAIR2::interpBadAntW(d, *i, near);
      const LibAIR2::InterpArrayData::wvrdata_t &data(d.g_wvrdata());
      for(size_t ii=0; ii<d.g_time().size(); ++ii){
	for(size_t k=0; k < 4; ++k){
	  double p=0;
	  for(LibAIR2::AntSetWeight::const_iterator j=near.begin(); j!=near.end(); ++j){
	    double thisData = data(ii,j->second,k);
	    if(thisData>0){
	      p+=thisData*j->first;
	    }
	    else{ // no good data; set solution to zero => will be flagged later
	      p=0.;
	      break;
	    }
	  }
	  d.set(ii, *i, k, p);
	}
      }
    }
    else
    { 
      *itsLog << LogIO::WARN << "Antenna " << *i 
	      << " has bad or no WVR and only " << near.size() << " near antennas ("
	      << maxdist_m << " m max. distance) to interpolate from. Required are " 
	      << minnumants << "." << LogIO::POST;

      for(size_t j=0; j<d.g_time().size(); ++j)
      {
        for(size_t k=0; k < 4; ++k)
        {
          d.set(j, *i, k, 0.); // set all WVR data of this antenna to zero => will be flagged later
        }
      }
      interpImpossibleAnts.insert(*i);
    }
  }
  
}


/// Work out which spectral windows might need to be reversed
std::set<size_t> reversedSPWs(const LibAIR2::MSSpec &sp,
			      const bool reverse,
			      const std::vector<long> &reversespw)
{
  std::set<size_t> to_be_reversed;
  if (reverse)
  {
    for (size_t spw =0; spw<sp.spws.size(); ++spw)
      to_be_reversed.insert(spw);
  }    
  if (reversespw.size() > 0)
  {
    std::vector<long> torev(reversespw);
    for(size_t i=0; i<torev.size(); ++i)
    {
      if (torev[i]<0 or torev[i] >= (int)sp.spws.size())
      {
	throw LibAIR2::SPWIDError(torev[i], 
				 sp.spws.size());
      }
      to_be_reversed.insert(torev[i]);
    }
  }    
  return to_be_reversed;
}

void printExpectedPerf(const LibAIR2::ArrayGains &g,
		       const LibAIR2::dTdLCoeffsBase &coeffs,
		       const std::vector<std::pair<double, double> > &tmask,
		       casacore::LogIO *itsLog)
{

  std::cout<<"  Expected performance "<<std::endl
	   <<"------------------------------------------------------------------"<<std::endl;
  
  std::vector<double> cr, err;
  coeffs.repr(cr, err);
  std::cout<<"* Estimated WVR thermal contribution to path fluctuations (micron per antenna): "
	   <<LibAIR2::thermal_error(cr)/1e-6
	   <<std::endl;
  const double grmsbl=g.greatestRMSBl(tmask);
  std::cout<<"* Greatest Estimated path fluctuation is (micron on a baseline): "
	   <<grmsbl/1e-6
	   <<std::endl;
  if (cr[1]>0.){
    std::cout<<"* Rough estimate path error due to coefficient error (micron on a baseline): "
	     <<grmsbl* (err[1]/cr[1])/1e-6
	     <<std::endl
	     <<std::endl;
  }
  else if(err[1]==0.) {
    std::cout<<"* Rough estimate path error due to coefficient error could not be calculated (is nominally zero)."<<std::endl;
  }
  else {
    std::cout<<"* Rough estimate path error due to coefficient error could not be calculated."<<std::endl;
    *itsLog << LogIO::WARN <<"* Rough estimate path error due to coefficient error could not be calculated." << LogIO::POST;
  }
}

/** Compute the time intervals over which the statistics should be
    computed
 */

void statTimeMask(const casacore::MeasurementSet &ms,
		  const std::vector<std::string> &statfield,
		  const std::vector<std::string> &statsource,
		  std::vector<std::pair<double, double> > &tmask,
		  const std::vector<size_t> &sortedI,
		  const std::vector<int> &wvrspws,
		  casacore::LogIO *itsLog)
{
  std::vector<int> flds;
  std::vector<double> time;
  std::vector<int> src;
  LibAIR2::fieldIDs(ms,
		   time,
		   flds,
		   src,
		   sortedI);
  std::vector<size_t> spws;
  LibAIR2::dataSPWs(ms, spws, sortedI);

  if (statfield.size() == 0 && statsource.size() == 0)
  {
    tmask.resize(0);
    tmask.push_back(std::pair<double, double>(time[0], time[time.size()-1]));
  }
  else if ( statsource.size() > 0)
  {
    std::set<size_t> fselect=LibAIR2::getSrcFields(ms,
						   statsource[0]);

    LibAIR2::fieldTimes(time,
		       flds,
		       spws,
		       fselect,
		       (size_t) wvrspws[0],
		       tmask);    
  }
  else
  {
    std::vector<std::string> fields(statfield);
    LibAIR2::field_t fnames=LibAIR2::getFieldNames(ms);

    std::set<size_t> fselect;
    size_t val;
    if (std::accumulate( fnames.begin( ),
			 fnames.end( ), false,
			 [&](bool acc, const LibAIR2::field_t::value_type &p) {
			   bool cmp = p.second == fields[0];
			   if ( cmp ) val = p.first;
			   return acc || cmp;
			 }
			 ))  // User supplied  field *name*
      {
	fselect.insert(val);
      }
    else
      {
	try 
	  {
	    size_t n=casacore::String::toInt(fields[0], true);
	    fselect.insert(n);
	  }
	catch (casacore::AipsError x)
	  {
	    *itsLog << LogIO::WARN <<"Could not understand statfield argument. Will use zeroth field."
		    << LogIO::POST;
	  }
      }
    LibAIR2::fieldTimes(time,
			flds,
			spws,
			fselect,
			(size_t) wvrspws[0],
			tmask);
  }
  LibAIR2::printStatTimes(std::cout,
			  time,
			  tmask);
}
		  

/// Compute the discrepance in path estimate between channels 1 and 3
void computePathDisc(const LibAIR2::InterpArrayData &d,
		     const std::vector<std::pair<double, double> > &tmask,
		     LibAIR2::dTdLCoeffsBase  &coeffs,
		     std::vector<double> &res)
{
  LibAIR2::ArrayGains g1(d.g_time(), 
			d.g_el(),
			d.g_state(),
			d.g_field(),
			d.g_source(),
			d.nAnts);
  std::array<double, 4> c1mask = {{0, 1, 0,0}};
  std::array<double, 4> c3mask = {{0, 0, 0,1}};
  std::array<double, 4> callmask = {{1, 1, 1,1}};
    
  coeffs.chmask=c1mask;
  g1.calc(d,
	  coeffs);    
  
  LibAIR2::ArrayGains g3(d.g_time(), 
			 d.g_el(),
			 d.g_state(),
			 d.g_field(),
			 d.g_source(),
			 d.nAnts);
  coeffs.chmask=c3mask;
  g3.calc(d,
	  coeffs);    

  g1.pathDiscAnt(g3, 
		 tmask,
		 res);

  coeffs.chmask=callmask;
  
}

void printFieldSegments(const std::vector<std::pair<double, double> >  &fb,
			double tbase)
{
  for (size_t i=0; i<fb.size(); ++i)
    {
      std::cout<<fb[i].first-tbase<<","<<fb[i].second-tbase<<std::endl;
    }

}

std::vector<std::set<std::string> > getTied(const std::vector<std::string> &input)
{

  std::vector<std::set<std::string> > res;

  for (size_t i=0; i< input.size(); ++i){
    std::set<std::string> cs;
    const std::string &par=input[i];
    char *pch;
    char *s = new char[par.size()+1];
    strcpy( s, par.c_str() );
    pch = strtok(s, ",");   
    while (pch != NULL){
      std::string first(pch);
      cs.insert(first);
      pch = strtok(NULL, ",");
    }
    delete [] s;
    res.push_back(cs);
  }

  return res;
}

// Convert tied source names to tied source IDs
std::vector<std::set<size_t> >  tiedIDs(const std::vector<std::set<std::string> > &tied,
					const casacore::MeasurementSet &ms,
					casacore::LogIO *itsLog)
{
  std::map<size_t, std::string > srcmap=LibAIR2::getSourceNames(ms);
  std::vector<std::set<size_t> > res;
  for (size_t i=0; i<tied.size(); ++i)
  {
    std::set<size_t> cs;
    for(std::set<std::string>::const_iterator j=tied[i].begin();
	j!=tied[i].end();
	++j)
    {	
      try
      {
	int srcid=casacore::String::toInt(*j, true);
	std::map<size_t, std::string>::const_iterator it = srcmap.find(srcid);
	if(it == srcmap.end()) { // id does not exist
	  *itsLog << LogIO::WARN << "Parameter 'tie': The source id " << *j << " is an integer but not a valid numerical Source ID. Will try to interpret it as a name ..." << LogIO::POST;
	  throw std::exception();
	}
	cs.insert(srcid);
      }
      catch (casacore::AipsError x)
      {
        size_t match;
        if ( std::accumulate( srcmap.begin( ),
                              srcmap.end( ), false,
                              [&](bool acc, const std::map<size_t, std::string >::value_type &p) {
				bool cmp = p.second == *j;
				if ( cmp ) match = p.first;
				return acc || cmp;
                              }
			      )) {
          cs.insert(match);
        } else {
          std::ostringstream oss;
          oss << "Parameter 'tie': The field " << *j << " is not recognised. Please check for typos." << std::endl;
          throw LibAIR2::WVRUserError(oss.str());
	}
      }
    } // end for
    res.push_back(cs);
  }
  return res;
}

void printTied(const std::vector<std::set<std::string> > &tied,
	       const std::vector<std::set<size_t> > &tiedi)
{
  for(size_t i=0; i<tied.size(); ++i)
  {
    std::set<std::string>::const_iterator it=tied[i].begin();
    std::cout << "Tying: " << *it;
    for(it++; it!=tied[i].end(); it++){
      std::cout <<" and "<<*it;
    }
    std::cout << std::endl;
  }
  if (tied.size())
    std::cout <<"Tied sets as numerical source IDs:"<<std::endl;
  for(size_t i=0; i<tiedi.size(); ++i)
  {
    std::set<size_t>::const_iterator it=tiedi[i].begin();
    std::cout << "Tying: " << *it;
    for(it++; it!=tiedi[i].end(); it++){
      std::cout <<" and "<<*it;
    }
    std::cout << std::endl;
  }

}

/** Compute the set of source_ids corresponding to a vector of source
    names.
 */
std::set<size_t> sourceSet(const std::vector<std::string> &sources,
			   const casacore::MeasurementSet &ms)
{
  std::map<size_t, std::string > snames=LibAIR2::getSourceNames(ms);
  std::set<size_t> sset;
  for(size_t i=0; i<sources.size(); ++i) {
	size_t match;
	if (std::accumulate( snames.begin( ),
			     snames.end( ), false,
			     [&](bool acc, const std::map<size_t, std::string >::value_type &p) {
			       bool cmp = p.second == sources[i];
			       if ( cmp ) match = p.first;
			       return acc || cmp;
			     }
			     )) {
	  sset.insert(match);
	}
  }
  return sset;
}
  

/** Filter the set of input WVR measurements to retrieve the
    coefficients from to exclude flagged sources

    This function takes and returns two containers: the first is the
    list of WVR values to be analysed; the second is a vector of time
    ranges to use.
    
    These time ranges must be filtered together with the inputs since
    they are not referenced directly to the inputs are simply assumed
    to be "row-synchronous".
    
 */
std::pair<LibAIR2::ALMAAbsInpL,  std::vector<std::pair<double, double> > >
filterInp(const LibAIR2::ALMAAbsInpL &inp,
	  const std::vector<std::pair<double, double> > &fb,
	  const std::vector<std::string> &sourceflag,
	  const casacore::MeasurementSet &ms)
{

  std::set<size_t> flagset=sourceSet(sourceflag, ms);

  LibAIR2::ALMAAbsInpL res;
  std::vector<std::pair<double, double> > rfb;
  size_t j=0;
  for(LibAIR2::ALMAAbsInpL::const_iterator i=inp.begin();
      i!=inp.end();
      ++i, ++j)
  {
    if (flagset.count(i->source) ==0)
    {
      res.push_back(*i);
      rfb.push_back(fb[j]);
    }
  }
  return std::make_pair(res, rfb);
}

/** Filter the set of input WVR measurements to retrieve the
    coefficients from to exclude flagged data points (zero Tobs)
    
 */
std::pair<LibAIR2::ALMAAbsInpL,  std::vector<std::pair<double, double> > >
filterFlaggedInp(const LibAIR2::ALMAAbsInpL &inp,
		 const std::vector<std::pair<double, double> > &fb)
{

  LibAIR2::ALMAAbsInpL res;
  std::vector<std::pair<double, double> > rfb;
  size_t j=0;
  bool fbFilled = (fb.size()>0);
  for(LibAIR2::ALMAAbsInpL::const_iterator i=inp.begin();
      i!=inp.end();
      ++i, ++j)
  {
    if(i->TObs[0]>0.) // flagged ALMAAbsInp would have TObs==0
    {
      res.push_back(*i);
      if(fbFilled)
      {
	rfb.push_back(fb[j]);
      }
    }
  }
  return std::make_pair(res, rfb);
}

/** Return the set of antenna IDs that do not have a WVR
 */
LibAIR2::AntSet NoWVRAnts(const LibAIR2::aname_t &an)
{
  LibAIR2::AntSet res;
  for(LibAIR2::aname_t::const_iterator i=an.begin();
      i!= an.end();
      ++i)
  {
    if (i->second[0]=='C' and i->second[1]=='M')
      res.insert(i->first);
  }
  return res;
}

std::string b_to_s(bool val){
  if(val){
    return("True");
  }
  else{
    return("False");
  }
}

void vl_to_os(const std::vector<long> &vl, std::ostringstream &os){
  if(vl.size()==0){
    os << "[]";
  }
  else{
    os << "[";
    std::copy(vl.begin(), vl.end(),
	      ostream_iterator<long>(os, ", "));
    os << "]";
  }
  return;
}    

void vs_to_os(const std::vector<std::string> &vs, std::ostringstream &os){
  if(vs.size()==0){
    os << "[]";
  }
  else{
    os << "[\"";
    std::copy(vs.begin(), vs.end(),
	      ostream_iterator<std::string>(os, "\", \""));
    os << "\"]";
  }
  return;
}    

namespace casac {

  /** default tool constructor 
   */
  wvr::wvr() {
    itsLog = new casacore::LogIO();
  }

  wvr::~wvr() {
    if(itsLog) delete itsLog;
    itsLog=0;
  }

  /** The WVRGCAL function
   */
  long wvr::gcal(const std::string &vis_par,
		 const std::string &output_par,
		 double toffset_par,
		 long nsol_par,
		 bool segsource_par,
		 bool reverse_par,
		 const std::vector<long> &reversespw_par,
		 bool disperse_par,
		 bool cont_par,
		 const std::string &wvrflag_par,
		 const std::vector<std::string> &sourceflag_par,
		 const std::vector<std::string> &statfield_par,
		 const std::vector<std::string> &statsource_par,
		 const std::vector<std::string> &tie_par,
		 long smooth_par,
		 double scale_par,
		 double maxdistm_par,
		 long minnumants_par,
		 double mingoodfrac_par,
		 bool usefieldtab_par,
		 const std::vector<long> &spw_par,
		 const std::vector<long> &wvrspw_par,
		 const std::string &refant_par,
		 const std::string &offsets_par,
		 const std::string &logfile_par
		 )
  {

    std::string logfile = logfile_par;
    if(logfile.length()<1){
      logfile = "wvrgcal.log";
    }
    
    casa::AsdmStMan::registerClass();

    int rval = -2;

    std::string cmdLineHistory; // to be entered into the gain table history
    {
      std::ostringstream toss;
      toss << "wvrgcal(vis=\""<<vis_par<<"\","<<
	" output=\""<<output_par<<"\","<<
	" toffset="<<std::to_string(toffset_par)<<","<<
	" nsol="<<to_string(nsol_par)<<","<<
	" segsource="<<b_to_s(segsource_par)<<","<<
	" reverse="<<b_to_s(reverse_par)<<","<<
	" reversespw=";
      vl_to_os(reversespw_par, toss);
      toss<<","<<
	" disperse="<<b_to_s(disperse_par)<<","<<
	" cont="<<b_to_s(cont_par)<<","<<
	" wvrflag=\""<<wvrflag_par<<"\","<<
	" sourceflag=";
      vs_to_os(sourceflag_par, toss);
      toss<<", statfield=";
      vs_to_os(statfield_par, toss);
      toss<<", statsource=";
      vs_to_os(statsource_par, toss);
      toss<<", tie=";
      vs_to_os(tie_par, toss);
      toss<<","<<
	" smooth="<<std::to_string(smooth_par)<<","<<
	" scale="<<std::to_string(scale_par)<<","<<
	" maxdistm="<<std::to_string(maxdistm_par)<<","<<
	" minnumants="<<std::to_string(minnumants_par)<<","<<
	" mingoodfrac="<<std::to_string(mingoodfrac_par)<<","<<
	" usefieldtab="<<b_to_s(usefieldtab_par)<<","<<
	" spw=[";
      vl_to_os(spw_par, toss);
      toss<<", wvrspw=";
      vl_to_os(wvrspw_par, toss);
      toss<<","<<
	" refant=\""<<refant_par<<"\","<<
	" offsets=\""<<offsets_par<<"\")";
      cmdLineHistory = toss.str();
    }
          
    // Check parameters

    if (vis_par.length()<1){
      *itsLog << LogIO::SEVERE << "No input measurement set given -- aborting " << LogIO::POST;
      return -1;
    }

    if (output_par.length()<1){
      *itsLog << LogIO::SEVERE << "No output file give -- aborting " << LogIO::POST;
      return -1;
    }

    if (segsource_par && cont_par){
      *itsLog << LogIO::SEVERE << "Multiple retrievals using continuum estimation not yet supported";
      *itsLog << "Set only one of  segfield and cont to True" << LogIO::POST;
      return -1;
    }

    if (nsol_par>1 && cont_par){
      *itsLog << LogIO::SEVERE << "Multiple retrievals using continuum estimation not yet supported";
      *itsLog << "You can not use the nsol parameter yet with cont" << LogIO::POST;
      return -1;
    }

    if (sourceflag_par.size()>0 && !segsource_par){
      *itsLog << LogIO::SEVERE << "Can only flag a source using sourceflag if the segsource option is also used";
      *itsLog << "Please either remove the sourceflag option or also specify the segsource option" << LogIO::POST;
      return -1;
    }

    if (tie_par.size()>0 && !segsource_par){
      *itsLog << LogIO::SEVERE << "Can only tie sources together if the segsource option is also used";
      *itsLog << "Please either remove the tie option or also specify the segsource option" << LogIO::POST;
      return -1;
    }

    if (reverse_par and reversespw_par.size()>0){
      *itsLog << LogIO::WARN << "You are specifying both the reverse and reversespw options. "
	"The latter will be ignored and all spectral windows will be reversed" << LogIO::POST;
    }

    if (smooth_par < 1){
      *itsLog << "smooth parameter must be 1 or greater" << LogIO::POST;
      return -1;
    }

    if (scale_par <= 0.){
      *itsLog << LogIO::SEVERE << "scale parameter must be > 0." << LogIO::POST;
      return -1;
    }

    if (maxdistm_par < 0.){
      *itsLog << LogIO::SEVERE << "maxdistm parameter must be 0. or greater" << LogIO::POST;
      return -1;
    }

    if (offsets_par.length() > 0){
      std::string offsetstable=offsets_par;
      try{
	casacore::Table f(offsetstable);
      }
      catch(const casacore::AipsError rE){
	*itsLog << LogIO::SEVERE << rE.getMesg() << std::endl;
	*itsLog << "The --offsets option needs to point to an existing CASA Table containing temperature offsets." << LogIO::POST;
	return -1;
      }
    }

    if (statfield_par.size() > 0){
      *itsLog << LogIO::WARN << "The use of \"statfield\" is not recommended as mosaiced"
	"observations are recorded as many separate fields. Recommended option"
	"to use is \"statsource\". " << LogIO::POST;
    }
  
    std::vector<std::set<std::string> > tied=getTied(tie_par);

    std::string msname(vis_par);
    casacore::MeasurementSet ms(msname);
    if (statsource_par.size()>0){
      std::string srcname=statsource_par[0];
      std::set<size_t> fselect=LibAIR2::getSrcFields(ms,
						     srcname);
      if (fselect.size() == 0){
	*itsLog<< LogIO::WARN << "No Source table entries appear to be identified with source "<< srcname <<std::endl
	       <<" that you supplied to the statsource option."<<std::endl
	       <<" Statistics will be corrupted and wvrgcal may fail."<< LogIO::POST;
      }
    }

    std::vector<int> wvrspws;
    {
      LibAIR2::SPWSet thewvrspws=LibAIR2::WVRSPWIDs(ms);
      if (wvrspw_par.size()>0){
	for(size_t i=0; i<wvrspw_par.size(); i++){
	  wvrspws.push_back((int) wvrspw_par[i]);
	  if(thewvrspws.count(wvrspws[i])==0){
	    *itsLog << LogIO::SEVERE << "SPW " << wvrspws[i] << " is not a WVR SPW or invalid." << LogIO::POST;
	    return -10;
	  }
	}
	*itsLog << LogIO::NORMAL << "Will use the following WVR SPWs:" <<std::endl;
	for(size_t i=0; i<wvrspws.size();i++){
	  *itsLog << " " << wvrspws[i];
	}
	*itsLog << LogIO::POST;
      }
      else{ // no wvrspws specified

	*itsLog << LogIO::NORMAL <<"Will use all WVR SPWs:"<<std::endl;
	for(LibAIR2::SPWSet::const_iterator si=thewvrspws.begin(); si!=thewvrspws.end();++si){
	  wvrspws.push_back((int) *si);
	  *itsLog << " " <<  *si;
	}
	*itsLog << LogIO::POST;
      }
    }

    std::vector<int> sciencespws;

    if (spw_par.size()>0){
      LibAIR2::SPWSet thewvrspws=LibAIR2::WVRSPWIDs(ms);

      int nspw=LibAIR2::numSPWs(ms);
      for(size_t i=0; i<spw_par.size(); i++){
	sciencespws.push_back((int) spw_par[i]);
	if(thewvrspws.count(sciencespws[i])!=0){
	  *itsLog << LogIO::WARN << "SPW "<< sciencespws[i] << " is a WVR SPW, not a science SPW." << LogIO::POST;
	}
	if(sciencespws[i]<0 || sciencespws[i]>= nspw){
	  *itsLog << LogIO::SEVERE <<"Invalid SPW "<< sciencespws[i] << LogIO::POST;
	  return -11;
	}
      }
      *itsLog << LogIO::NORMAL << "Will produce solutions for the following SPWs:"<<std::endl;
      for(size_t i=0; i<sciencespws.size();i++){
	*itsLog << " " << sciencespws[i];
      }
      *itsLog << LogIO::POST;
    }
    else{ // no sciencespws specified
      *itsLog << LogIO::NORMAL <<"Will produce solutions for all SPWs:"<<std::endl;
      for(size_t i=0; i<LibAIR2::numSPWs(ms);i++){
	sciencespws.push_back(i);
	*itsLog << " " << i;
      }
      *itsLog << LogIO::POST;
    }

    std::string fnameout=output_par;

    std::string offsetstable=offsets_par;

    std::set<size_t> useID=LibAIR2::skyStateIDs(ms);

    LibAIR2::AntSet wvrflagset;
    // Prepare flagging and interpolation
    if (wvrflag_par.length()>0){
      wvrflagset=getAntPars(wvrflag_par, ms);    
    }

    LibAIR2::aname_t anames=LibAIR2::getAName(ms);
    LibAIR2::AntSet nowvr=NoWVRAnts(anames);
  
    LibAIR2::AntSet interpwvrs(wvrflagset); // the antennas to interpolate solutions for
    interpwvrs.insert(nowvr.begin(), nowvr.end());

    LibAIR2::AntSet flaggedants; // the antennas flagged in the ANTENNA table are not to be interpolated
    LibAIR2::WVRAddFlaggedAnts(ms, flaggedants);

    wvrflagset.insert(flaggedants.begin(),flaggedants.end());

    if(interpwvrs.size()+flaggedants.size()==ms.antenna().nrow()){
      *itsLog << LogIO::SEVERE << "No good antennas with WVR data found." << LogIO::POST;
      return -1;
    }

    // redirect stdout to logfile
    std::fstream logstream; 
    logstream.open(logfile, std::ios::out);
    std::streambuf* stream_buffer_cout_orig = std::cout.rdbuf();
    std::streambuf* stream_buffer_logstream = logstream.rdbuf();
    std::cout.rdbuf(stream_buffer_logstream);

    LibAIR2::printBanner(std::cout);
    
    int iterations = 0;

    try{

      while(rval<0 && iterations<2){

	iterations++;

	std::vector<size_t> sortedI; // to be filled with the time-sorted row number index
	std::set<int> flaggedantsInMain; // the antennas totally flagged in the MS main table
	std::unique_ptr<LibAIR2::InterpArrayData> d (LibAIR2::loadWVRData(ms,
									  wvrspws,
									  sortedI, 
									  flaggedantsInMain,
									  mingoodfrac_par,
									  usefieldtab_par==false, // i.e. usepointing==true
									  offsetstable)
						     );
      
	// For debug purposes, print the loaded WVR data: 
	// for(size_t j=0; j<d->g_time().size(); ++j)
	// {
	// 	for(size_t i=0; i < ms.antenna().nrow(); ++i)
	// 	  {
	// 	    std::cout << "row ant data " << j << " " << i << " "
	// 		      << d->g_wvrdata()[j][i][0] << " "
	// 		      << d->g_wvrdata()[j][i][1] << " " 
	// 		      << d->g_wvrdata()[j][i][2] << " "
	// 		      << d->g_wvrdata()[j][i][3] << std::endl; 
	// 	  }
	// }
     

	interpwvrs.insert(flaggedantsInMain.begin(),flaggedantsInMain.end()); // for flagInterp()
	wvrflagset.insert(flaggedantsInMain.begin(),flaggedantsInMain.end());

	d->offsetTime(toffset_par);
     
	if (smooth_par > 1){
	  smoothWVR(*d, (int) smooth_par);
	}
     
	d.reset(LibAIR2::filterState(*d, useID));

	LibAIR2::AntSet interpImpossibleAnts;

	// Flag and interpolate
	flagInterp(ms,
		   interpwvrs,
		   *d,
		   maxdistm_par,
		   (int) minnumants_par,
		   interpImpossibleAnts,
		   itsLog);

	// Determine the reference antenna for dTdL calculation
	int therefant = -1; 

	if (refant_par.length()>0){
	  std::vector<size_t> refants=getAntParsV(refant_par, ms);    
	  for(std::vector<size_t>::iterator it=refants.begin(); it != refants.end(); it++){ // 
	    if(interpImpossibleAnts.count(*it)==0){
	      therefant = *it; // use the first of the given list of possible ref antennas which was OK or which could be interpolated to
	      break;
	    }
	    else{
	      std::cout << "Given reference antenna " << *it << "==" << anames.at(*it) 
			<< " is flagged and cannot be interpolated." << std::endl;
	    }	   
	  }
	  if(therefant<0){
	    std::cout << "None of the given reference antennas is usable." << std::endl;
	    *itsLog << LogIO::WARN << "None of the given reference antennas is usable." << LogIO::POST;
	    rval = -1;
	    break; // leave while loop
	  }

	}
	else{
	  LibAIR2::AntSet wvrants=LibAIR2::WVRAntennas(ms, wvrspws);
	  for(LibAIR2::AntSet::iterator it=wvrants.begin(); it != wvrants.end(); it++){
	    if(interpImpossibleAnts.count(*it)==0){
	      therefant = *it; // use the first antenna which was OK or which could be interpolated to
	      break;
	    }
	  }
	  if(therefant<0){
	    std::cout << "No antennas with sufficient WVR data found." << std::endl;
	    *itsLog << LogIO::WARN << "No antennas with sufficient WVR data found." << std::endl;
	    rval = -1;
	    break; // leave while loop
	  }
	}

	std::cout << "Choosing";
	if(interpwvrs.count(therefant)>0){
	  std::cout << " (interpolated)";
	}
	std::cout << " antenna " << therefant  << " == " << anames.at(therefant)
		  << " as reference antenna for dTdL calculations." << std::endl;


	LibAIR2::ArrayGains g(d->g_time(), 
			      d->g_el(),
			      d->g_state(),
			      d->g_field(),
			      d->g_source(),
			      d->nAnts);
     
	std::unique_ptr<LibAIR2::dTdLCoeffsBase>  coeffs;
     
	// These are the segments on which coefficients are re-calculated
	std::vector<std::pair<double, double> >  fb;
     
	if ( cont_par ){
	  std::cout<<"[Output from \"cont\" option has not yet been updated]"
		   <<std::endl;
	  coeffs.reset(LibAIR2::SimpleSingleCont(*d, therefant));
	}
	else{
	
	  LibAIR2::ALMAAbsInpL inp;
	  if(segsource_par){
	    std::vector<int> flds;
	    std::vector<double> time;
	    std::vector<int> src;
	    LibAIR2::fieldIDs(ms, 
			      time,
			      flds,
			      src,
			      sortedI);
	    try{
	      std::vector<std::set<size_t> >  tiedi=tiedIDs(tied, ms, itsLog);
	      
	      printTied(tied, tiedi);
	      LibAIR2::fieldSegmentsTied(time,
					 src,
					 tiedi,
					 fb);
	    }
	    catch(LibAIR2::WVRUserError& x){
	      std::cout << x.what() << std::endl;
	      *itsLog << LogIO::WARN << x.what() << LogIO::POST;
	      rval = -1;
	      break; // leave while loop
	    }

	    //printFieldSegments(fb, time[0]);
	    
	    //       { // debugging output
	    // 	std::vector<double> tt(d->g_time());
	    // 	std::vector<double> te(d->g_el());
	    // 	std::vector<size_t> ts(d->g_state());
	    // 	std::vector<size_t> tf(d->g_field());
	    // 	std::vector<size_t> tsou(d->g_source());
	    // 	for(uint i=0; i<tt.size(); i++){
	    // 	  *itsLog << LogIO::WARN << "i time el state field source " << i << " " << tt[i] << " ";
	    // 	  *itsLog << LogIO::WARN << te[i] << " ";
	    // 	  *itsLog << LogIO::WARN << ts[i] << " ";
	    // 	  *itsLog << LogIO::WARN << tf[i] << " ";
	    // 	  *itsLog << LogIO::WARN << tsou[i] << std::endl;
	    // 	}
	    // 	*itsLog << LogIO::WARN << "nAnts " << d->nAnts << std::endl;
	    // 	for(uint i=0; i<time.size(); i++){
	    // 	  *itsLog << LogIO::WARN << "i time  " << i << " " << time[i] << std::endl;
	    // 	}
	    // 	for(uint i=0; i<fb.size(); i++){
	    // 	  *itsLog << LogIO::WARN << "i fb  " << i << " " << fb[i].first << " " << fb[i].second  << std::endl;
	    // 	}
	    // 	for(std::set<size_t>::iterator it = useID.begin(); it != useID.end(); it++){
	    // 	  *itsLog << LogIO::WARN << "useID " << *it << std::endl;
	    // 	}
	    //       }
	    
	    inp=FieldMidPointI(*d,
			       fb,
			       useID,
			       therefant);
	    
	  }
	  else{
	    const size_t n=nsol_par;
	    inp=LibAIR2::MultipleUniformI(*d, 
					  n,
					  useID,
					  therefant);
	  }

	
	  if (sourceflag_par.size() > 0){
	    std::tie(inp,fb)=filterInp(inp,
				       fb,
				       sourceflag_par,
				       ms);
	  }
	
	  std::tie(inp,fb)=filterFlaggedInp(inp,
					    fb);

	  *itsLog << LogIO::NORMAL << "Calculating the coefficients now ... " << LogIO::POST;
	  std::cerr << "Calculating the coefficients now ... " << std::endl;
	  LibAIR2::ALMAResBaseList rlist;
	  LibAIR2::AntSet problemAnts;

	  rval = 0;

	  try {
	    rlist=LibAIR2::doALMAAbsRet(inp,
					fb,
					problemAnts);
	  }
	  catch(const std::runtime_error& rE){
	    rval = 1;
	    *itsLog << LogIO::WARN << "problem while calculating coefficients:"
		    << std::endl << "         LibAIR2::doALMAAbsRet: " << rE.what() << LogIO::POST;
	    std::cout << std::endl << "WARNING: problem while calculating coefficients:"
		      << std::endl << "         LibAIR2::doALMAAbsRet: " << rE.what() << std::endl;
	  }
	
	  if(problemAnts.size()>0){
	    
	    rval = -2;
	    
	    if(iterations<2){
	      for(LibAIR2::AntSet::const_iterator it=problemAnts.begin(); it!=problemAnts.end(); it++){
		if(interpwvrs.count(*it)==0){
		  *itsLog << LogIO::WARN << "Flagging antenna " << *it << " == " << anames.at(*it) << LogIO::POST;
		  std::cout	<< "Flagging antenna " << *it << " == " << anames.at(*it) << std::endl;
		  interpwvrs.insert(*it); // for flagInterp()
		  wvrflagset.insert(*it); // for later log output
		}
	      }
	      *itsLog << LogIO::WARN << "Reiterating ..." << LogIO::POST;
	      std::cout	<< "Reiterating ..." << std::endl;
	      continue;
	    }
	    else{
	      *itsLog << LogIO::WARN << "Number of remaining antennas with problematic WVR measurements: " << problemAnts.size() << LogIO::POST;
	      std::cout << "Number of remaining antennas with problematic WVR measurements: " << problemAnts.size() << std::endl;
	      *itsLog << LogIO::WARN << "Will continue without further iterations ..." << LogIO::POST;
	      std::cout << "Will continue without further iterations ..." << std::endl;
	    }	      
	  }	   
	
	  std::cerr<<"done!"
		   <<std::endl;
	
	  std::cout<<"       Retrieved parameters      "<<std::endl
		   <<"----------------------------------------------------------------"<<std::endl
		   << rlist.ptr_list <<std::endl;
	  
	  if (segsource_par){
	    //std::vector<int> flds;
	    //std::vector<double> time;
	    //std::vector<int> src;
	    //LibAIR2::fieldIDs(ms, 
	    //		    time,
	    //		    flds,
	    //		    src,
	    //		    sortedI);
	    
	    coeffs.reset(LibAIR2::SimpleMultiple(fb,
						   rlist));   
	  }
	  else{
	    coeffs.reset(LibAIR2::ALMAAbsProcessor(inp, rlist));
	  }  
	  
	}
    
	try{
	  g.calc(*d,
		 *coeffs);    
	}
	catch(const std::runtime_error& rE){
	  std::cout << "Problem while calculating gains: " << rE.what() << std::endl;
	  *itsLog << LogIO::WARN << "Problem while calculating gains: " << rE.what() << LogIO::POST;
	  rval = 1;
	  break; // leave while loop
	}

	if (sourceflag_par.size() > 0){
	  std::set<size_t> flagset=sourceSet(sourceflag_par,
					     ms);
	  g.blankSources(flagset);
	}
     
	std::vector<std::pair<double, double> > tmask;
	statTimeMask(ms, statfield_par, statsource_par, tmask, sortedI, wvrspws, itsLog);
	
	std::vector<double> pathRMS;
	g.pathRMSAnt(tmask, pathRMS);
	
	std::vector<double> pathDisc;
	try{
	  computePathDisc(*d, 
			  tmask,
			  *coeffs,
			  pathDisc);
     
	  std::cout<<LibAIR2::AntITable(anames,
					wvrflagset,
					nowvr,
					pathRMS,
					pathDisc,
					interpImpossibleAnts);
       
	  printExpectedPerf(g, 
			    *coeffs,
			    tmask,
			    itsLog);
       
	}
	catch(const std::runtime_error& rE){
	  std::cout << "Problem while calculating path RMS discrepancy: " << rE.what() << std::endl;
	  *itsLog << LogIO::WARN << "Problem while calculating path RMS discrepancy: " << rE.what() << LogIO::POST;
	  rval = 1;
	  break; // leave while loop
	}
     
	if (scale_par != 1.0){
	  g.scale(scale_par);
	}
     
	LibAIR2::MSSpec sp;
	loadSpec(ms, sciencespws, sp);
	std::set<size_t> to_be_reversed=reversedSPWs(sp, reverse_par, reversespw_par);  
     
	std::cout << "Writing gain table ..." << std::endl;

	// Write new table, including history
	LibAIR2::writeNewGainTbl(g,
				 fnameout.c_str(),
				 sp,
				 to_be_reversed,
				 disperse_par,
				 msname,
				 cmdLineHistory,
				 interpImpossibleAnts);


      } // end while
    }
    catch(const std::runtime_error& rE){
      std::cout << "Problem while processing WVR data: " << rE.what() << std::endl;
      *itsLog << LogIO::WARN << "Problem while processing WVR data: " << rE.what() << LogIO::POST;
      rval = 1;
    }
    
    std::cout.rdbuf(stream_buffer_cout_orig);
    logstream.close();
    return rval;
  } // end gcal
}
