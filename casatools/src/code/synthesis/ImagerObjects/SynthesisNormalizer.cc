//# SynthesisNormalizer.cc: Implementation of Gather/Scatter for parallel major cycles
//# Copyright (C) 1997-2008
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This program is free software; you can redistribute it and/or modify it
//# under the terms of the GNU General Public License as published by the Free
//# Software Foundation; either version 2 of the License, or (at your option)
//# any later version.
//#
//# This program is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
//# more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with this program; if not, write to the Free Software Foundation, Inc.,
//# 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: casa-feedback@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
//# $Id$

#include <casacore/casa/Exceptions/Error.h>
#include <iostream>
#include <sstream>

#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>

#include <casacore/casa/Logging.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Logging/LogMessage.h>
#include <casacore/casa/Logging/LogSink.h>
#include <casacore/casa/Logging/LogMessage.h>

#include <casacore/casa/OS/DirectoryIterator.h>
#include <casacore/casa/OS/File.h>
#include <casacore/casa/OS/Path.h>

#include <casacore/casa/OS/HostInfo.h>
#include <casacore/lattices/Lattices/LatticeLocker.h>
#include <casacore/images/Images/TempImage.h>
#include <casacore/images/Images/SubImage.h>
#include <casacore/images/Regions/ImageRegion.h>

#include <synthesis/ImagerObjects/SynthesisNormalizer.h>

#include <sys/types.h>
#include <unistd.h>
using namespace std;

using namespace casacore;
namespace casa { //# NAMESPACE CASA - BEGIN
  
  SynthesisNormalizer::SynthesisNormalizer() :
      itsImages(std::shared_ptr<SIImageStore>()),
      itsPartImages(Vector<std::shared_ptr<SIImageStore> >()),
      itsImageName(""),
      itsPartImageNames(Vector<String>(0)),
      itsPBLimit(0.2),
      itsMapperType("default"),
      itsNTaylorTerms(1),
      itsNFacets(1),
      itsIsSingleDish(False)
  {
    itsFacetImageStores.resize(0);
    itsPBLimit = 0.35;
  }
  
  SynthesisNormalizer::~SynthesisNormalizer() 
  {
    LogIO os( LogOrigin("SynthesisNormalizer","destructor",WHERE) );
    os << LogIO::DEBUG1 << "SynthesisNormalizer destroyed" << LogIO::POST;
    SynthesisUtilMethods::getResource("End SynthesisNormalizer");

  }
  
  
  void SynthesisNormalizer::setupNormalizer(Record normpars)
  {
    LogIO os( LogOrigin("SynthesisNormalizer","setupNormalizer",WHERE) );

    try {
      if ( normpars.isDefined("psfcutoff") ) { // A single string
        normpars.get( RecordFieldId("psfcutoff") , itsPsfcutoff );
      } else {
        itsPsfcutoff=0.35;
      }

      if ( normpars.isDefined("imagename") ) { // A single string
        itsImageName = normpars.asString( RecordFieldId("imagename")); 
      } else {
        throw AipsError("imagename not specified");
      }

      if ( normpars.isDefined("partimagenames") ) { // A vector of strings
        normpars.get( RecordFieldId("partimagenames") , itsPartImageNames );
      } else {
        itsPartImageNames.resize(0);
      }

      if ( normpars.isDefined("pblimit") ) {
        normpars.get( RecordFieldId("pblimit") , itsPBLimit );
      } else {
        itsPBLimit = 0.2;
      }

      if ( normpars.isDefined("normtype") ) { // A single string
        itsNormType = normpars.asString( RecordFieldId("normtype"));
      } else {
        itsNormType = "flatnoise"; // flatnoise, flatsky
      }

      { // Debug comments
        //      cout << "Chosen normtype : " << itsNormType << endl;

        // For multi-term choices.
        // Try to eliminate, after making imstores hold aux descriptive info.
        /*
        if ( normpars.isDefined("mtype") )  // A single string
          { itsMapperType = normpars.asString( RecordFieldId("mtype")); }
        else
          { itsMapperType = "default";}
        */
      }

      if ( normpars.isDefined("deconvolver") ) { // A single string
        String dec = normpars.asString( RecordFieldId("deconvolver"));
        if (dec == "mtmfs") { itsMapperType="multiterm"; }
        else itsMapperType="default";
      } else {
        itsMapperType = "default";
      }

      if ( normpars.isDefined("nterms") ) { // A single int
        itsNTaylorTerms = normpars.asuInt( RecordFieldId("nterms")); }
      else {
        itsNTaylorTerms = 1;
      }

      if ( normpars.isDefined("facets") ) { // A single int
        itsNFacets = normpars.asuInt( RecordFieldId("facets"));
      } else { 
        itsNFacets = 1;
      }

      if ( normpars.isDefined("restoringbeam") ) {
        if (normpars.dataType("restoringbeam") == TpString) {
          itsUseBeam = normpars.asString( RecordFieldId("restoringbeam") );
        } else {
          itsUseBeam = "";
        } 
      } // else: keep itsUseBeam unchanged

      if ( normpars.isDefined("makesingledishnormalizer") ) { // A single bool
        itsIsSingleDish =
          normpars.asBool( RecordFieldId("makesingledishnormalizer") );
      } else { 
        itsIsSingleDish = False;
      }

    }
    catch (AipsError &x) {
      throw AipsError(
          "Error in reading gather/scatter parameters: " + x.getMesg()
        );
    }

  } // end of setupNormalizer


  void SynthesisNormalizer::gatherImages(Bool dopsf, Bool doresidual, Bool dodensity)
  {
    //    cout << " partimagenames :" << itsPartImageNames << endl;

    Bool needToGatherImages = setupImagesOnDisk();
    //if(dopsf)
    //  cerr << "GATHER " << itsImages->psf()->shape() << endl;

    if( needToGatherImages )
      {
	LogIO os( LogOrigin("SynthesisNormalizer", "gatherImages",WHERE) );

	AlwaysAssert( itsPartImages.nelements()>0 , AipsError );
	//Bool doresidual = !dopsf;
        Bool doweight = dopsf ; //|| ( doresidual && ! itsImages->hasSensitivity() );
        //	Bool doweight = dopsf || ( doresidual && itsImages->getUseWeightImage(*(itsPartImages[0]->residual())) );
	
	os << "Gather "<< (doresidual?"residual":"") << ( (dopsf&&doresidual)?",":"")  
	   << (dopsf?"psf":"") << ( (dopsf&&doweight)?",":"")  
	   << (doweight?"weight":"")<< " images : " << itsPartImageNames 
	   << " onto :" << itsImageName << LogIO::POST;
	
	// Add intelligence to modify all only the first time, but later, only residual;
	itsImages->resetImages( /*psf*/dopsf, /*residual*/doresidual, /*weight*/doweight ); 
	
	for( uInt part=0;part<itsPartImages.nelements();part++)
	  {
	    itsImages->addImages( itsPartImages[part], /*psf*/dopsf, /*residual*/doresidual, /*weight*/doweight, /*griddedwt*/dodensity );
	    itsPartImages[part]->releaseLocks();
	  }

      }// end of image gathering.
    
    // Normalize by the weight image.
    //    divideResidualByWeight();
    itsImages->releaseLocks();
    
  }// end of gatherImages

  void SynthesisNormalizer::gatherPB()
  {
    if( itsPartImageNames.nelements()>0 )
      {

	try
	  {
	    itsPartImages[0] = makeImageStore( itsPartImageNames[0] );
	  }
	catch(AipsError &x)
	  {
	    throw(AipsError("Cannot construct ImageStore for "+itsPartImageNames[0]+". "+x.what()));
	  }

	try{
	    LatticeExpr<Float> thepb( *(itsPartImages[0]->pb()) );
	    LatticeLocker lock1(*(itsImages->pb()), FileLocker::Write);
	    itsImages->pb()->copyData(thepb);

	  }
	catch(AipsError &x)
	  {
	    throw(AipsError("Cannot copy the PB : "+x.getMesg()));
	  }

      }
  }

  void SynthesisNormalizer::scatterModel()
  {

    LogIO os( LogOrigin("SynthesisNormalizer", "scatterModel",WHERE) );

    setupImagesOnDisk(); // To open up and initialize itsPartImages.

    //    os << "In ScatterModel : " << itsPartImages.nelements() << " for " << itsPartImageNames << LogIO::POST;
    
    if( itsPartImages.nelements() > 0 ) // && itsImages->doesImageExist(modelName) )
      {
	os << "Send the model from : " << itsImageName << " to all nodes :" << itsPartImageNames << LogIO::POST;
	
	// Make the list of model images. This list is of length >1 only for multi-term runs.
	// Vector<String> modelNames( itsImages->getNTaylorTerms() );
	// if( modelNames.nelements() ==1 ) modelNames[0] = itsImages->getName()+".model";
	// if( modelNames.nelements() > 1 ) 
	//   {
	//     for( uInt nt=0;nt<itsImages->getNTaylorTerms();nt++)
	//       modelNames[nt] = itsImages->getName()+".model.tt" + String::toString(nt);
	//   }


	Vector<String> modelNames( itsImages->getNTaylorTerms() );
	if( itsImages->getType()=="default" ) modelNames[0] = itsImages->getName()+".model";
	if( itsImages->getType()=="multiterm" ) 
	  {
	    for( uInt nt=0;nt<itsImages->getNTaylorTerms();nt++)
	      modelNames[nt] = itsImages->getName()+".model.tt" + String::toString(nt);
	  }
	
	
	
	for( uInt part=0;part<itsPartImages.nelements();part++)
	  {
	    itsPartImages[part]->setModelImage( modelNames );
	    itsPartImages[part]->releaseLocks();
	  }
	itsImages->releaseLocks();
      }
    
  }// end of scatterModel
  
//----------------

void SynthesisNormalizer::gatherWeightDensity(){

  Bool needToGatherImages = setupImagesOnDisk();
  // if(dopsf)
  //   cerr << "GATHER " << itsImages->psf()->shape() << endl;

  if (needToGatherImages) {
    LogIO os(LogOrigin("SynthesisNormalizer", "gatherWeightDensity", WHERE));

    itsImages->gridwt()->set(0.0);
    Record iminfo(itsImages->gridwt()->miscInfo());
    Vector<Float> d2out;
    Vector<Float> f2out;
    if(iminfo.isDefined("d2")){
      iminfo.get("d2", d2out);
      iminfo.get("f2", f2out);
    }
    Vector<Float> f2part(f2out.nelements(), 0.0);
    bool recalcF2 = (d2out.nelements() >0 && d2out[0] == 1.0); // need to recalculate F2 for normalized Briggs
    if(recalcF2)
      f2out.set(0.0);
    Int nxim = itsImages->gridwt()->shape()(0);
    Int nyim = itsImages->gridwt()->shape()(1);
    for (uInt part = 0; part < itsPartImages.nelements(); part++) {
      itsImages->addImages(itsPartImages[part], false,
                           /*residual*/ false, /*weight*/ false,
                           /*griddedwt*/ true);
      Record partiminfo(itsPartImages[part]->gridwt()->miscInfo());
      mergeWeightDensityInfo(iminfo, partiminfo);
      if (recalcF2) {
        Vector<Float> f2in;
        (itsPartImages[part]->gridwt()->miscInfo()).get("f2", f2in);
        Array<Float> arr;
        if (d2out.nelements() > 1) {
          IPosition blc(5, 0, 0, 0, 0, 0);
          IPosition shp(5, nxim, nyim, 1, 1, 1);
          for (uInt k = 0; k < d2out.nelements(); ++k) {
            blc[4] = k;
            itsPartImages[part]->gridwt()->getSlice(arr, blc, shp, True);
            f2in[k] = f2in[k] * sum(arr * arr);
          }
        }
        else{ //not multifield weighting
          itsPartImages[part]->gridwt()->get(arr, True);
          f2in[0] = f2in[0] * sum(arr * arr);
        }
        f2out += f2in;
      }
      itsPartImages[part]->releaseLocks();
    }
    if(recalcF2){
      Array<Float> arr;
      if (d2out.nelements() > 1) {
        IPosition blc(5, 0, 0, 0, 0, 0);
        IPosition shp(5, nxim, nyim, 1, 1, 1);
        for (uInt k = 0; k < d2out.nelements(); ++k) {
          blc[4] = k;
          itsImages->gridwt()->getSlice(arr, blc, shp, True);
          Float sumarr2 = sum(arr * arr);
          if(sumarr2> 0)
            f2out[k] = f2out[k] / sumarr2;
        }
      }
      else{
        itsImages->gridwt()->get(arr, True);
        Float sumarr2 = sum(arr * arr);
        if (sumarr2 > 0)
          f2out[0] = f2out[0] / sumarr2;
      }
      iminfo.define("f2", f2out);
      
    }
    itsImages->gridwt()->setMiscInfo(iminfo);
    
  } // end of image gathering.

  // Normalize by the weight image.
  //    divideResidualByWeight();
  itsImages->releaseLocks();

}//end of gatherweightdensity

//-------------------------


  string SynthesisNormalizer::scatterWeightDensity()
  {

    LogIO os( LogOrigin("SynthesisNormalizer", "scatterWeightDensity",WHERE) );

    setupImagesOnDisk(); // To open up and initialize itsPartImages.
    string weightname = "";
    //    os << "In ScatterModel : " << itsPartImages.nelements() << " for " <<
    //    itsPartImageNames << LogIO::POST;

    if( itsPartImages.nelements() > 0 )
      {
	os << "Send the gridded weight from : " << itsImageName << " to all nodes :" << itsPartImageNames << LogIO::POST;
	
	for( uInt part=0;part<itsPartImages.nelements();part++)
	  {
	    itsPartImages[part]->setWeightDensity( itsImages );
	    itsPartImages[part]->releaseLocks();
	  }
          weightname=itsImages->gridwt()->name();
          itsImages->releaseLocks();
      }
      return weightname;
  } // end of gatherImages

  void SynthesisNormalizer::divideResidualByWeight()
  {
    LogIO os( LogOrigin("SynthesisNormalizer", "divideResidualByWeight",WHERE) );

    if( itsNFacets==1) {
      itsImages->divideResidualByWeight( itsPBLimit, itsNormType );
    }
    else {
      for ( uInt facet=0; facet<itsNFacets*itsNFacets; facet++ )
        { itsFacetImageStores[facet]->divideResidualByWeight( itsPBLimit , itsNormType ); }
    }
    itsImages->releaseLocks();
    
  }
  
  void SynthesisNormalizer::divideResidualByWeightSD()
  {
    LogIO os( LogOrigin("SynthesisNormalizer", "divideResidualByWeightSD",WHERE) );

    if( itsNFacets==1) {
      itsImages->divideResidualByWeightSD( itsPBLimit );
    }
    else {
      for ( uInt facet=0; facet<itsNFacets*itsNFacets; facet++ )
        { itsFacetImageStores[facet]->divideResidualByWeightSD( itsPBLimit ); }
    }
    itsImages->releaseLocks();

  }
  void SynthesisNormalizer::makePSFBeamset(){
    LogIO os(LogOrigin("SynthesisNormalizer", "dividePSFByWeight", WHERE));
    try{
      if(!itsImages){
        itsImages = makeImageStore( itsImageName, false );
      }


    }
    catch(AipsError& x){

      throw(AipsError("Programmers error no psf is made on disk and  trying to fit"));
    }
    itsImages->makeImageBeamSet(itsPsfcutoff);
    itsImages->releaseLocks();
  }


  void SynthesisNormalizer::divideWeightBySumWt() {

    LogIO os(LogOrigin("SynthesisNormalizer", "divideWeightBySumWt", WHERE));
    itsImages->divideWeightBySumwt();

  }
  

  
  void SynthesisNormalizer::dividePSFByWeight()
  {
    LogIO os( LogOrigin("SynthesisNormalizer", "dividePSFByWeight",WHERE) );
    {
     
      LatticeLocker lock1 (*(itsImages->psf()), FileLocker::Write);
    
      if( itsNFacets==1) {
        itsImages->dividePSFByWeight(itsPBLimit);
      }
      else {
        // Since PSFs are normed just by their max, this sequence is OK.
        setPsfFromOneFacet();
        itsImages->dividePSFByWeight(itsPBLimit);
      }
    }// release of lock1 otherwise releaselock below will cause it to core
    //dump as the psf pointer goes dangling
      // Check PSF quality by fitting beams
    {
      itsImages->calcSensitivity();
    
      itsImages->makeImageBeamSet(itsPsfcutoff);
      Bool verbose(False);
      if (itsUseBeam=="common") verbose=True;
      itsImages->printBeamSet(verbose);
    }
    itsImages->releaseLocks();

  }

  void SynthesisNormalizer::normalizePrimaryBeam()
  {
    LogIO os( LogOrigin("SynthesisNormalizer", "normalizePrimaryBeam",WHERE) );

    if( itsImages==NULL ) { itsImages = makeImageStore( itsImageName ); }

    gatherPB();

    // Irrespective of facets.
    itsImages->normalizePrimaryBeam(itsPBLimit);
  }


  void SynthesisNormalizer::divideModelByWeight()
  {
    LogIO os( LogOrigin("SynthesisNormalizer", "divideModelByWeight",WHERE) );
    if( ! itsImages ) 
      {
	//os << LogIO::WARN << "No imagestore yet. Trying to construct, so that the starting model can be divided by wt if needed...." << LogIO::POST;
	
	try
	  {
	    itsImages = makeImageStore( itsImageName );
	  }
	catch(AipsError &x)
	  {
	    throw(AipsError("Cannot construct ImageStore for "+itsImageName+". "+x.what()));
	  }
	//	return;
      }
    if( itsNFacets==1) {
      itsImages->divideModelByWeight( itsPBLimit, itsNormType );
    }
    else {
      for ( uInt facet=0; facet<itsNFacets*itsNFacets; facet++ )
        { itsFacetImageStores[facet]->divideModelByWeight( itsPBLimit, itsNormType ); }
    }
 }

  void SynthesisNormalizer::multiplyModelByWeight()
  {
    //    LogIO os( LogOrigin("SynthesisNormalizer", "multiplyModelByWeight",WHERE) );
    if( itsNFacets==1) {
      itsImages->multiplyModelByWeight( itsPBLimit , itsNormType );
    }
    else {
      for ( uInt facet=0; facet<itsNFacets*itsNFacets; facet++ )
        { itsFacetImageStores[facet]->multiplyModelByWeight( itsPBLimit , itsNormType); }
    }
 }


  std::shared_ptr<SIImageStore> SynthesisNormalizer::getImageStore()
  {
    LogIO os( LogOrigin("SynthesisNormalizer", "getImageStore", WHERE) );
    return itsImages;
  }

  void SynthesisNormalizer::setImageStore( SIImageStore* imstore )
  {
    LogIO os( LogOrigin("SynthesisNormalizer", "setImageStore", WHERE) );
    itsImages.reset( imstore );
  }

  void SynthesisNormalizer::setImageStore( std::shared_ptr<SIImageStore>& imstore )
  {
    LogIO os( LogOrigin("SynthesisNormalizer", "setImageStore", WHERE) );
    itsImages= imstore ;
    
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////    Internal Functions start here.  These are not visible to the tool layer.
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  Bool SynthesisNormalizer::setupImagesOnDisk() {
    LogIO os( LogOrigin("SynthesisNormalizer","setupImagesOnDisk",WHERE) );

    Bool needToGatherImages = false;

    String err("");

    // Check if full images exist, and open them if possible.
    Bool foundFullImage = false;
    try {
      itsImages = makeImageStore( itsImageName );
      foundFullImage = true;
    }
    catch (AipsError &x) {
      //throw( AipsError("Error in constructing a Deconvolver : "+x.getMesg()) );
      err = err += String(x.getMesg()) + "\n";
      foundFullImage = false;
    }

    os << LogIO::DEBUG2
      << " Found full images : " << foundFullImage
      << LogIO::POST;

    // Check if part images exist
    Bool foundPartImages = itsPartImageNames.nelements() > 0 ? true : false ;
    itsPartImages.resize( itsPartImageNames.nelements() );

    for ( uInt part=0; part < itsPartImageNames.nelements() ; part++ ) {
      try {
        itsPartImages[part] = makeImageStore ( itsPartImageNames[part] );
        foundPartImages |= true;
      }
      catch (AipsError &x) {
        // throw AipsError(
        //        "Error in constructing a Deconvolver : " + x.getMesg()
        //       );
        err = err += String(x.getMesg()) + "\n";
        foundPartImages = false;
      }
    }

    os << LogIO::DEBUG2
      << " Found part images : " << foundPartImages
      << LogIO::POST;

    if ( not foundPartImages ) {
      if ( foundFullImage and itsPartImageNames.nelements() > 0 ) {
        // Pick the coordsys, etc from fullImage,
        // and construct new/fresh partial images.
        os << LogIO::DEBUG2
          << "Found full image, but no partial images."
          << " Making partImStores for : " << itsPartImageNames
          << LogIO::POST;

        String imopen = itsImages->getName() + ".residual"
          + ( (itsMapperType == "multiterm") ? ".tt0" : "" );
        Directory imdir( imopen );
        if ( not imdir.exists() ) {
          imopen = itsImages->getName() + ".psf"
          + ( (itsMapperType == "multiterm") ? ".tt0" : "" );
          Directory imdir2( imopen );
          if ( not imdir2.exists() )
            throw AipsError(
              "Cannot find partial image psf or residual for "
              + itsImages->getName() + err
            );
        }

        PagedImage<Float> temppart(imopen);

        Bool useweightimage =
          itsImages->getUseWeightImage( *(itsImages->sumwt()) );
        for ( uInt part=0; part<itsPartImageNames.nelements(); part++ ) {
          itsPartImages[part] = makeImageStore(
            itsPartImageNames[part], temppart, useweightimage);
        }
        foundPartImages = True;
      } else {
        itsPartImages.resize(0); 
        foundPartImages = False;
      }
    } else { // Check that all have the same shape.
      AlwaysAssert( itsPartImages.nelements() > 0 , AipsError );
      IPosition tempshape = itsPartImages[0]->getShape();
      for ( uInt part=1; part<itsPartImages.nelements(); part++ ) {
        if ( tempshape != itsPartImages[part]->getShape() ) {
          throw AipsError(
            "Shapes of partial images to be combined, do not match" + err
          );
        }
      }
    }

    // Make sure all images exist and are consistent with each other.
    // At the end, itsImages should be valid
    if ( foundPartImages ) { // Partial Images exist. Check that 'full' exists,
                             // and do the gather. 
      if ( foundFullImage ) { // Full image exists.
                              // Just check that shapes match with parts.
        os << LogIO::DEBUG2
          << "Partial and Full images exist."
          << "Checking if part images have the same shape as the full image : ";

        IPosition fullshape = itsImages->getShape();
        for ( uInt part=0; part < itsPartImages.nelements() ; part++ ) {
          IPosition partshape = itsPartImages[part]->getShape();
          if ( partshape != fullshape ) {
            os << LogIO::DEBUG2
              << "NO" << LogIO::POST;
            throw AipsError(
              "Shapes of the partial and full images on disk do not match."
              " Cannot gather" + err
            );
          }
        }
        os << LogIO::DEBUG2 << "Yes" << LogIO::POST;

      } else { // Full image does not exist. 
               // Need to make it, using the shape and coords of part[0]
        os << LogIO::DEBUG2 
          << "Only partial images exist. Need to make full images"
            << LogIO::POST  ;

        AlwaysAssert( itsPartImages.nelements() > 0, AipsError );

        // Find an image to open and pick csys,shape from.
        String imopen = itsPartImageNames[0] + ".residual"
          + ( (itsMapperType == "multiterm") ? ".tt0" : "" );
        Directory imdir( imopen );
        if ( not imdir.exists() ) {
          imopen = itsPartImageNames[0] + ".psf"
            + ( (itsMapperType == "multiterm") ? ".tt0" : "" );
          Directory imdir2( imopen );
          if ( not imdir2.exists() ) {
            imopen = itsPartImageNames[0] + ".gridwt";
            Directory imdir3( imopen );
            if ( not imdir3.exists() )
              throw AipsError(
                "Cannot find partial image psf or residual or gridwt for "
                + itsPartImageNames[0] + err
              );
          }
        }

        PagedImage<Float> temppart( imopen );

        Bool useweightimage =
          itsPartImages[0]->getUseWeightImage( *(itsPartImages[0]->sumwt()) );
        itsImages = makeImageStore (itsImageName, temppart, useweightimage);
        foundFullImage = true;
      }

      // By now, all partial images and the full images exist on disk,
      // and have the same shape.
      needToGatherImages = true;
    } else { // No partial images supplied. Operating only with full images.

      if ( foundFullImage ) {
          os << LogIO::DEBUG2
            << "Full images exist : " << itsImageName
            << LogIO::POST;
      }
      else { // No full image on disk either. Error.
          throw AipsError(
            "No images named " + itsImageName 
            + " found on disk. No partial images found either." + err
          );
      }
    }

    // Remove ?
    itsImages->psf();

    itsImages->validate();

    // Set up facet Imstores..... if needed
    if ( itsNFacets > 1 ) {
      // First, make sure that full images have been allocated
      // before trying to make references.....
      //    if ( not itsImages->checkValidity(
      //                true/*psf*/,
      //                true/*res*/,
      //                true/*wgt*/,
      //                true/*model*/,
      //                false/*image*/,
      //                false/*mask*/,
      //                true/*sumwt*/ ) ) {
      //      throw AipsError(
      //        "Internal Error : Invalid ImageStore for "
      //         + itsImages->getName())
      //      );
      //    }

      //        Array<Float> ttt;
      //        (itsImages->sumwt())->get(ttt);
      //        cout << "SUMWT full : " << ttt <<  endl;
        
      itsFacetImageStores.resize( itsNFacets*itsNFacets );
      for ( uInt facet=0; facet<itsNFacets*itsNFacets; ++facet ) {
        itsFacetImageStores[facet] =
          itsImages->getSubImageStore(facet, itsNFacets);
        // Array<Float> qqq;
        // itsFacetImageStores[facet]->sumwt()->get(qqq);
        // cout << "SUMWTs for " << facet << " : " << qqq << endl;
      }
    }

    os << LogIO::DEBUG2
      << "Need to Gather ? " << needToGatherImages
      << LogIO::POST;

    itsImages->releaseLocks();

    return needToGatherImages;
  } // end of setupImagesOnDisk


  std::shared_ptr<SIImageStore> SynthesisNormalizer::makeImageStore(
    const String &imagename,
    const bool useweightimage)
  {
    // The constructors use ignoresumwt so use the !useweightimage
    if ( itsMapperType == "multiterm" ) {
      return std::shared_ptr<SIImageStore>(
        new SIImageStoreMultiTerm(
          imagename,
          itsNTaylorTerms,
          /* ignorefacets= */ true,
          /* ignoresumwt=  */ not useweightimage
        )
      );
    } else {
      return std::shared_ptr<SIImageStore>(
        new SIImageStore(
          imagename,
          /* ignorefacets= */ true,
          /* noRequireSumwt= */ not useweightimage,
          /* makeSingleDishStore= */ itsIsSingleDish
        )
      );
    }
    itsImages->releaseLocks();
  }


 /**
  * build a new ImageStore, whether SIImageStore or SIImageStoreMultiTerm, borrowing
  * image information from one partial image.
  *
  * @param imagename image name for the new SIStorageManager
  * @param part partial image from which miscinfo, etc. data will be borrowed
  * @param useweightimage useweight option for the new SIStorageManager
  *
  * @return A new SIImageStore object for the image name given.
  */
  std::shared_ptr<SIImageStore> SynthesisNormalizer::makeImageStore(const String &imagename,
                                                               const PagedImage<Float> &part,
                                                               Bool useweightimage)
  {
    // borrow shape, coord, imageinfo and miscinfo
    auto shape = part.shape();
    auto csys = part.coordinates();
    auto objectname = part.imageInfo().objectName();
    auto miscinfo = part.miscInfo();
    if (itsMapperType == "multiterm") {
      std::shared_ptr<SIImageStore> multiTermStore =
          std::make_shared<SIImageStoreMultiTerm>(
              imagename, csys, shape, objectname, miscinfo, itsNFacets, false,
              itsNTaylorTerms, useweightimage);
      return multiTermStore;
    } else {
      return std::make_shared<SIImageStore>(imagename, csys, shape, objectname,
                                            miscinfo, false, useweightimage);
    }
  }

  //
  //---------------------------------------------------------------
  //
   void SynthesisNormalizer::setPsfFromOneFacet()
   {
       // Copy the PSF from one facet to the center of the full image, to use for the minor cycle
    //
    // This code segment will work for single and multi-term 
    // - for single term, the index will always be 0, and SIImageStore's access functions know this.
    // - for multi-term, the index will be the taylor index and SIImageStoreMultiTerm knows this.

     /// itsImages
     /// itsFacetImageStores

      {
	IPosition shape=(itsImages->psf(0))->shape();
	IPosition blc(4, 0, 0, 0, 0);
	IPosition trc=shape-1;
	for(uInt tix=0; tix<2 * itsImages->getNTaylorTerms() - 1; tix++)
	  {
	    TempImage<Float> onepsf((itsFacetImageStores[0]->psf(tix))->shape(), 
				    (itsFacetImageStores[0]->psf(tix))->coordinates());
	    onepsf.copyData(*(itsFacetImageStores[0]->psf(tix)));
	    //now set the original to 0 as we have a copy of one facet psf
	    (itsImages->psf(tix))->set(0.0);
	    blc[0]=(shape[0]-(onepsf.shape()[0]))/2;
	    trc[0]=onepsf.shape()[0]+blc[0]-1;
	    blc[1]=(shape[1]-(onepsf.shape()[1]))/2;
	    trc[1]=onepsf.shape()[1]+blc[1]-1;
	    Slicer sl(blc, trc, Slicer::endIsLast);
	    SubImage<Float> sub(*(itsImages->psf(tix)), sl, true);
	    sub.copyData(onepsf);
	  }
      }

      //      cout << "In setPsfFromOneFacet : sumwt : " << itsImages->sumwt()->get() << endl;
  
   }

   void SynthesisNormalizer::mergeWeightDensityInfo(Record& outinfo, const Record& partinfo){
     uInt mapsize=0;
     if(outinfo.isDefined("multimapsize"))
      outinfo.get("multimapsize", mapsize);
     uInt partsize;
     if(!partinfo.isDefined("multimapsize"))
       return;
     partinfo.get("multimapsize", partsize);
     std::vector<String> key2add;
     std::vector<Int> val2add;
     for (uInt k = 0; k < partsize; ++k){
      String key;
      partinfo.get("key" + String::toString(k), key);
      Bool hasKey = False;
      for (uInt j = 0; j < mapsize; ++j) {
        String imkey;
        outinfo.get("key" + String::toString(j), imkey);
        hasKey = hasKey || (imkey == key);
      }  
      if (!hasKey) {
          key2add.push_back(key);
          Int val;
          partinfo.get("val" + String::toString(k), val);
          val2add.push_back(val);
        }
       }
     if (key2add.size() > 0) {
       for (uInt k = mapsize; k < (mapsize + key2add.size()); ++k) {
         String key = "key" + String::toString(k);
        outinfo.define(key, key2add[k - mapsize]);
        String val = "val" + String::toString(k);
        outinfo.define(val, val2add[k - mapsize]);
      }
      mapsize += key2add.size();
      outinfo.define("multimapsize", mapsize);
     }
   }

} //# NAMESPACE CASA - END

