//# table_convert.h
//# Copyright (C) 2022,2023
//# Associated Universities, Inc. Washington DC, USA.
//#
//# This library is free software; you can redistribute it and/or modify it
//# under the terms of the GNU Library General Public License as published by
//# the Free Software Foundation; either version 2 of the License, or (at your
//# option) any later version.
//#
//# This library is distributed in the hope that it will be useful, but WITHOUT
//# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
//# License for more details.
//#
//# You should have received a copy of the GNU Library General Public License
//# along with this library; if not, write to the Free Software Foundation,
//# Inc., 675 Massachusetts Ave, Cambridge, MA 02139, USA.
//#
//# Correspondence concerning AIPS++ should be addressed as follows:
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
#include <map>
#include <cstring>
#include <functional>
#include <Python.h>
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <casacore/casa/BasicSL/Complex.h>
#include <casacore/casa/BasicSL/String.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Containers/ValueHolder.h>

// ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
// These conversion functions are used by tablerow_cmpt.cc and table_getcoliter.cc to
// convert casacore C++ values into python values. These functions allocate python
// memory so the PyGILState_Ensure( ) must be used before calling these functions to
// ensure that the GIL is locked (preventing any other Python activity) while these
// functions do the conversion, allocating python ojects in the process.
//
// Previously, these functions were included as static functions within
// tablerow_cmpt.cc which was introduced as part of CAS-13894.
// 
// ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
namespace casac {

    inline size_t non_zero( size_t val ) {
        return  val <= 0 ? 5 : val;
    }

    // convert a boolean value to a PyObject
    inline PyObject *toPy( bool b ) {
        if ( b ) { Py_INCREF(Py_True); return Py_True; }
        else { Py_INCREF(Py_False); return Py_False; }
    }

    // convert numeric scalars to a PyObject
#define PY_NUM_SCALAR( CASACORE_TYPE, NUMPY_TYPE )                                 \
    inline PyObject *toPy( CASACORE_TYPE i ) {                                     \
        static PyObject *itemlen = PyLong_FromLong(sizeof(i));                     \
        return PyArray_Scalar( &i, PyArray_DescrFromType(NUMPY_TYPE), itemlen );   \
    }

    PY_NUM_SCALAR( int8_t, NPY_INT8 )
    PY_NUM_SCALAR( uint8_t, NPY_UINT8 )
    PY_NUM_SCALAR( int16_t, NPY_INT16 )
    PY_NUM_SCALAR( uint16_t, NPY_UINT16 )
    PY_NUM_SCALAR( int32_t, NPY_INT32 )
    PY_NUM_SCALAR( uint32_t, NPY_UINT32 )
    PY_NUM_SCALAR( int64_t, NPY_INT64 )
    PY_NUM_SCALAR( uint64_t, NPY_UINT64 )
    PY_NUM_SCALAR( float, NPY_FLOAT )
    PY_NUM_SCALAR( double, NPY_DOUBLE )
    PY_NUM_SCALAR( casacore::Complex, NPY_COMPLEX64 )
    PY_NUM_SCALAR( casacore::DComplex, NPY_COMPLEX128 )

    // convert a string to a PyObject
    inline PyObject *toPy( const casacore::String &s ) { return PyUnicode_FromString(s.c_str( )); }

    // convert an array of strings to a PyObject
    inline PyObject *toPy( const casacore::Array<casacore::String> &a ) {
        auto shape = a.shape( );
        size_t stringlen = std::accumulate( a.begin( ), a.end( ), (size_t) 0, []( size_t tally, const casacore::String &s ) { return s.size( ) > tally ? s.length( ) : tally; } );
        size_t memlen = a.nelements( ) * non_zero(stringlen) * sizeof(uint32_t);
        void *mem = PyDataMem_NEW(memlen);
        uint32_t *ptr = reinterpret_cast<uint32_t*>(mem);
        for ( const auto &str : a ) {
            for ( size_t i=0; i < non_zero(stringlen); ++i ) {
                *ptr++ = i < str.size( ) ? (unsigned char) str[i] : 0;
            }
        }
        return PyArray_New( &PyArray_Type, shape.nelements( ), (npy_intp*) shape.storage( ), NPY_UNICODE, nullptr, mem, non_zero(stringlen)*sizeof(uint32_t), NPY_ARRAY_OWNDATA | NPY_ARRAY_FARRAY, nullptr );
    }

    // convert numeric arrays to PyObjects
    // --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    // Allocating the result object with:
    //
    //   PyObject *ndarray = PyArray_New( &PyArray_Type, shape.nelements( ), (npy_intp*) shape.storage( ), NUMPY_TYPE,
    //                                    nullptr, nullptr, 0, NPY_ARRAY_FARRAY, nullptr );
    //
    // and then filling it after the fact with:
    //
    //   bool free_storage = false;
    //   auto storage = a.getStorage( free_storage );
    //   std::memcpy( PyArray_DATA( reinterpret_cast<PyArrayObject *>(ndarray)),
    //                              storage, a.nelements( ) * sizeof(CASACORE_TYPE) );
    //   PyArray_ENABLEFLAGS( reinterpret_cast<PyArrayObject *>(ndarray), NPY_ARRAY_OWNDATA );
    //   if ( free_storage ) delete storage;
    //
    // worked on RHEL7 + python 3.6 but fails on macos 10.15 + python 3.8
    //
#define PY_NUM_ARRAY( CASACORE_TYPE, NUMPY_TYPE )                                                                         \
    inline PyObject *toPy( const casacore::Array<CASACORE_TYPE> &a ) {                                                    \
        auto shape = a.shape( );                                                                                          \
        size_t memlen = a.nelements( ) * sizeof(CASACORE_TYPE);                                                           \
        auto *mem = PyDataMem_NEW(memlen);                                                                                \
        auto *ptr = reinterpret_cast<CASACORE_TYPE*>(mem);                                                                \
        for ( const auto &ele : a ) {                                                                                     \
            *ptr++ = ele;                                                                                                 \
        }                                                                                                                 \
        auto result = PyArray_New( &PyArray_Type, shape.nelements( ), (npy_intp*) shape.storage( ), NUMPY_TYPE, nullptr,  \
                                   mem, sizeof(CASACORE_TYPE), NPY_ARRAY_OWNDATA | NPY_ARRAY_FARRAY, nullptr );           \
        /*** setting NPY_ARRAY_OWNDATA here is required to avoid memory leak of allocated data (rhel7 + python 3.6) ***/  \
        PyArray_ENABLEFLAGS( reinterpret_cast<PyArrayObject*>(result), NPY_ARRAY_OWNDATA );                               \
        return result;                                                                                                    \
    }

    PY_NUM_ARRAY( bool, NPY_BOOL )
    PY_NUM_ARRAY( int8_t, NPY_INT8 )
    PY_NUM_ARRAY( uint8_t, NPY_UINT8 )
    PY_NUM_ARRAY( int16_t, NPY_INT16 )
    PY_NUM_ARRAY( uint16_t, NPY_UINT16 )
    PY_NUM_ARRAY( int32_t, NPY_INT32 )
    PY_NUM_ARRAY( uint32_t, NPY_UINT32 )
    PY_NUM_ARRAY( int64_t, NPY_INT64 )
    PY_NUM_ARRAY( uint64_t, NPY_UINT64 )
    PY_NUM_ARRAY( float, NPY_FLOAT )
    PY_NUM_ARRAY( double, NPY_DOUBLE )
    PY_NUM_ARRAY( casacore::Complex, NPY_COMPLEX64 )
    PY_NUM_ARRAY( casacore::DComplex, NPY_COMPLEX128 )

    inline PyObject *toPy( const casacore::Record &rec ) {
        using namespace casacore;

        // build map from table cell types to conversion functions
        std::map<int,std::function<PyObject*(size_t i)>> function_map = { {TpBool,[&](size_t i) ->PyObject* { return toPy(rec.asBool(i)); }},
                                                                          {TpChar,[&](size_t i) ->PyObject* { return toPy(rec.asuChar(i)); }},
                                                                          {TpUChar,[&](size_t i) ->PyObject* { return toPy(rec.asuChar(i)); }},
                                                                          {TpShort,[&](size_t i) ->PyObject* { return toPy(rec.asShort(i)); }},
                                                                          {TpUShort,[&](size_t i) ->PyObject* { return toPy(rec.asShort(i)); }},
                                                                          {TpInt,[&](size_t i) ->PyObject* { return toPy(rec.asInt(i)); }},
                                                                          {TpUInt,[&](size_t i) ->PyObject* { return toPy(rec.asuInt(i)); }},
                                                                          {TpInt64,[&](size_t i) ->PyObject* { return toPy((int64_t)rec.asInt64(i)); }},
                                                                          {TpFloat,[&](size_t i) ->PyObject* { return toPy(rec.asFloat(i)); }},
                                                                          {TpDouble,[&](size_t i) ->PyObject* { return toPy(rec.asDouble(i)); }},
                                                                          {TpComplex,[&](size_t i) ->PyObject* { return toPy(rec.asComplex(i)); }},
                                                                          {TpDComplex,[&](size_t i) ->PyObject* { return toPy(rec.asDComplex(i)); }},
                                                                          {TpArrayBool,[&](size_t i) ->PyObject* { return toPy(rec.asArrayBool(i)); }},
                                                                          {TpArrayUChar,[&](size_t i) ->PyObject* { return toPy(rec.asArrayuChar(i)); }},
                                                                          {TpArrayChar,[&](size_t i) ->PyObject* { return toPy(rec.asArrayuChar(i)); }},
                                                                          {TpArrayShort,[&](size_t i) ->PyObject* { return toPy(rec.asArrayShort(i)); }},
                                                                          {TpArrayUShort,[&](size_t i) ->PyObject* { return toPy(rec.asArrayShort(i)); }},
                                                                          {TpArrayInt,[&](size_t i) ->PyObject* { return toPy(rec.asArrayInt(i)); }},
                                                                          {TpArrayUInt,[&](size_t i) ->PyObject* { return toPy(rec.asArrayuInt(i)); }},
                                                                          {TpArrayFloat,[&](size_t i) ->PyObject* { return toPy(rec.asArrayFloat(i)); }},
                                                                          {TpArrayDouble,[&](size_t i) ->PyObject* { return toPy(rec.asArrayDouble(i)); }},
                                                                          {TpArrayComplex,[&](size_t i) ->PyObject* { return toPy(rec.asArrayComplex(i)); }},
                                                                          {TpArrayDComplex,[&](size_t i) ->PyObject* { return toPy(rec.asArrayDComplex(i)); }},
                                                                          {TpString,[&](size_t i) ->PyObject* { return toPy(rec.asString(i)); }},
                                                                          {TpRecord,[&](size_t i) ->PyObject* { return toPy(rec.asRecord(i)); }}
                                                                        };

        // create result
        auto result = PyDict_New( );
        if ( result == nullptr ) throw PyExc_MemoryError;
        // loop through record fields
        for ( uInt i=0; i < rec.nfields( ); ++i ) {
            auto func = function_map.find( rec.dataType(i) );
            // lookup conversion function
            if ( func != function_map.end( ) ) {
                auto newobj = func->second(i);
                auto name = PyUnicode_FromString(rec.name(i).c_str( ));
                // set field in result
                if ( PyDict_SetItem( result, name, newobj ) != 0 ) {
                    Py_DECREF(result);
                    Py_DECREF(newobj);
                    Py_DECREF(name);
                    throw PyExc_ValueError;
                }
                Py_DECREF(newobj);
                Py_DECREF(name);
            } else {
                Py_DECREF(result);
                throw PyExc_TypeError;
            }
        }
        return result;
    }

    inline PyObject *toPy( const casacore::ValueHolder &hldr ) {
        using namespace casacore;

        // build map from table cell types to conversion functions
        std::map<int,std::function<PyObject*( )>> function_map = { {TpBool,[&]( ) ->PyObject* { return toPy(hldr.asBool( )); }},
                                                                   {TpChar,[&]( ) ->PyObject* { return toPy(hldr.asuChar( )); }},
                                                                   {TpUChar,[&]( ) ->PyObject* { return toPy(hldr.asuChar( )); }},
                                                                   {TpShort,[&]( ) ->PyObject* { return toPy(hldr.asShort( )); }},
                                                                   {TpUShort,[&]( ) ->PyObject* { return toPy(hldr.asuShort( )); }},
                                                                   {TpInt,[&]( ) ->PyObject* { return toPy(hldr.asInt( )); }},
                                                                   {TpUInt,[&]( ) ->PyObject* { return toPy(hldr.asuInt( )); }},
                                                                   {TpInt64,[&]( ) ->PyObject* { return toPy((int64_t)hldr.asInt64( )); }},
                                                                   {TpFloat,[&]( ) ->PyObject* { return toPy(hldr.asFloat( )); }},
                                                                   {TpDouble,[&]( ) ->PyObject* { return toPy(hldr.asDouble( )); }},
                                                                   {TpComplex,[&]( ) ->PyObject* { return toPy(hldr.asComplex( )); }},
                                                                   {TpDComplex,[&]( ) ->PyObject* { return toPy(hldr.asDComplex( )); }},
                                                                   {TpString,[&]( ) ->PyObject* { return toPy(hldr.asString( )); }},
                                                                   {TpArrayBool,[&]( ) ->PyObject* { return toPy(hldr.asArrayBool( )); }},
                                                                   {TpArrayUChar,[&]( ) ->PyObject* { return toPy(hldr.asArrayuChar( )); }},
                                                                   {TpArrayChar,[&]( ) ->PyObject* { return toPy(hldr.asArrayuChar( )); }},
                                                                   {TpArrayShort,[&]( ) ->PyObject* { return toPy(hldr.asArrayShort( )); }},
                                                                   {TpArrayUShort,[&]( ) ->PyObject* { return toPy(hldr.asArrayuShort( )); }},
                                                                   {TpArrayInt,[&]( ) ->PyObject* { return toPy(hldr.asArrayInt( )); }},
                                                                   {TpArrayUInt,[&]( ) ->PyObject* { return toPy(hldr.asArrayuInt( )); }},
                                                                   {TpArrayInt64,[&]( ) ->PyObject* { return toPy(hldr.asArrayInt64( )); }},
                                                                   {TpArrayFloat,[&]( ) ->PyObject* { return toPy(hldr.asArrayFloat( )); }},
                                                                   {TpArrayDouble,[&]( ) ->PyObject* { return toPy(hldr.asArrayDouble( )); }},
                                                                   {TpArrayComplex,[&]( ) ->PyObject* { return toPy(hldr.asArrayComplex( )); }},
                                                                   {TpArrayDComplex,[&]( ) ->PyObject* { return toPy(hldr.asArrayDComplex( )); }},
                                                                   {TpArrayString,[&]( ) ->PyObject* { return toPy(hldr.asArrayString( )); }},
                                                                   {TpRecord,[&]( ) ->PyObject* { return toPy(hldr.asRecord( )); }}
                                                                 };
        
        auto func = function_map.find( hldr.dataType( ) );
        // lookup conversion function
        if ( func != function_map.end( ) ) {
            return func->second( );
        } else {
            throw PyExc_TypeError;
        }
    }
}
