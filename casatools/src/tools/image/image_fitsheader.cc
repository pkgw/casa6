#include <string.h>
#include <image_cmpt.h>
#include <casacore/images/Images/ImageFITSConverter.h>

#define _ORIGIN casacore::LogOrigin(_class, __func__, WHERE)

namespace casac {

    inline std::complex<double> cvtcomplex( const casacore::IComplex &c ) {
        return std::complex<double>( c.real( ), c.imag( ) );
    }
    inline std::complex<double> cvtcomplex( const casacore::Complex &c ) {
        return std::complex<double>( c.real( ), c.imag( ) );
    }
    inline std::complex<double> cvtcomplex( const casacore::DComplex &c ) {
        return std::complex<double>( c.real( ), c.imag( ) );
    }

    inline record kw2record( FitsKeywordList &kwl, const std::vector<std::string> &exclude ) {
        record result;
        for ( casacore::FitsKeyword *kv=(kwl.first( ),kwl.next( ));
              kv; kv=kwl.next( ) ) {
            if ( std::find(exclude.begin( ),exclude.end( ),kv->name( )) == std::end(exclude) ) {
                auto existing = result.find( kv->name() );
                switch (kv->type( )) {
                    case casacore::FITS::LOGICAL:
                        if ( existing == std::end(result) )
                            result.insert( kv->name(), variant(kv->asBool( )) );
                        else
                            existing->second.push(kv->asBool( ));
                        break;
                    case casacore::FITS::BYTE:
                    case casacore::FITS::SHORT:
                        // kv->asInt( ) --->> only, ONLY works for FITS::LONG
                        if ( existing == std::end(result) )
                            result.insert( kv->name(), variant((long)kv->asFloat( )) );
                        else
                            existing->second.push((long)kv->asFloat( ));
                        break;
                    case casacore::FITS::LONG:
                        if ( existing == std::end(result) )
                            result.insert( kv->name(), variant((long)kv->asInt( )) );
                        else
                            existing->second.push((long)kv->asInt( ));
                        break;
                    case casacore::FITS::FLOAT:
                        if ( existing == std::end(result) )                    
                            result.insert( kv->name(), variant(kv->asFloat( )) );
                        else
                            existing->second.push(kv->asFloat( ));
                        break;
                    case casacore::FITS::DOUBLE:
                        if ( existing == std::end(result) )
                            result.insert( kv->name(), variant(kv->asDouble( )) );
                        else
                            existing->second.push(kv->asDouble( ));
                        break;
                    case casacore::FITS::COMPLEX:
                        if ( existing == std::end(result) )
                            result.insert( kv->name(), variant(cvtcomplex(kv->asComplex( ))) );
                        else
                            existing->second.push(kv->asComplex( ));
                        break;
                    case casacore::FITS::ICOMPLEX:
                        if ( existing == std::end(result) )
                            result.insert( kv->name(), variant(cvtcomplex(kv->asIComplex( ))) );
                        else
                            existing->second.push(cvtcomplex(kv->asIComplex( )));
                        break;
                    case casacore::FITS::DCOMPLEX:
                        if ( existing == std::end(result) )
                            result.insert( kv->name(), variant(cvtcomplex(kv->asDComplex( ))) );
                        else
                            existing->second.push(cvtcomplex(kv->asDComplex( )));
                        break;
                    case casacore::FITS::STRING:
                        if ( existing == std::end(result) )
                            result.insert( kv->name(), variant(kv->asString( )) );
                        else
                            existing->second.push(kv->asString( ));
                        break;
                    case casacore::FITS::NOVALUE:
                        {
                            auto len = kv->commlen( ) <= 72 ? kv->commlen( ) : 72;
                            char buffer[73];
                            memcpy( buffer, kv->comm( ), len );
                            buffer[len] = '\0';

                            if ( existing == std::end(result) )
                                result.insert( kv->name(), variant(buffer) );
                            else
                                existing->second.push(std::string(buffer));
                        }
                        break;
                    case casacore::FITS::BIT:
                    case casacore::FITS::CHAR:
                    case casacore::FITS::VADESC:
                    case casacore::FITS::FSTRING:
                    case casacore::FITS::REAL:
                        if ( existing == std::end(result) ) {
                            record rec;
                            // no way to fetch any of these...
                            rec.insert("ERROR", std::string("cannot retrieve type ") + 
                                       ( kv->type( ) == casacore::FITS::BIT ? "BIT" :
                                         kv->type( ) == casacore::FITS::CHAR ? "CHAR" :
                                         kv->type( ) == casacore::FITS::VADESC ? "VADESC" :
                                         kv->type( ) == casacore::FITS::FSTRING ? "FSTRING" :
                                         kv->type( ) == casacore::FITS::REAL ? "REAL" : "<unknown-type>" ) );
                            result.insert( kv->name(), variant(rec) );
                        }
                        break;
                    default:
                        result.insert( kv->name(), variant(record( )) );
                }
            }
        }
        return result;
    }

    variant *image::fitsheader( bool retstr, const std::vector<std::string> &exclude ) {

        _log << _ORIGIN;
        if (_detached()) {
            return retstr ? new variant(std::string( )) :
                            new variant(record( ));
        }
        try {
            if ( ! _imageF ) {
                ThrowCc("fitsheader only supports floating point images");   
            }
            casacore::String error;
            casacore::ImageFITSHeaderInfo fhi;
            auto success = ImageFITSConverter::ImageHeaderToFITS( error, fhi, *_imageF );
            if ( success == false ) {
                ThrowCc("header creation failed: " + error);
            }
            size_t name_length = _imageF->name( ).length( );
            if ( name_length > 0 ) {
                //------------------------------------------------------
                //-- FITS header values are limited to 68 characters  --
                //------------------------------------------------------
                char *imagenme = strdup(_imageF->name( ).c_str( ));
                fhi.kw.mk( "IMAGENME", strlen(imagenme) > 68 ? &imagenme[name_length-68] : imagenme );
                free(imagenme);
            }

            // exclude not currently supported with retstr == true
            return retstr ? new variant(fhi.kw.toString( )) :
                new variant(kw2record(fhi.kw,exclude));
        } catch (const AipsError& x) {
            _log << LogIO::SEVERE << "Exception Reported: " << x.getMesg()
                 << LogIO::POST;
            RETHROW(x);
        }
        return retstr ? new variant(std::string( )) :
                        new variant(record( ));
        
    }
}


    
