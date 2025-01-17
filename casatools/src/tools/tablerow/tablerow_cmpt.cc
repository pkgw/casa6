//# tablerow_cmpt.cc
//# Copyright (C) 2022
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
//#        Internet email: casa-feedback@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#

#include <swigconvert_python.h>
#include <table_convert.h>
#include <table_cmpt.h>

#include <stdint.h>
#include <iostream>
#include <casacore/tables/Tables/TableProxy.h>
#include <casacore/tables/Tables/TableRowProxy.h>
#include <stdcasa/StdCasa/CasacSupport.h>
#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/Containers/Record.h>
#include <tablerow_cmpt.h>

using std::function;
using casa::toRecord;
using casa::fromRecord;
using casacore::LogIO;
using casacore::LogOrigin;
using casacore::AipsError;
using casacore::TableRowProxy;
using casacore::Array;
using casacore::Record;
using casacore::Complex;
using casacore::DComplex;
using casacore::String;

static bool _tablerow_initialize_numpy( ) {
    static bool initialized = false;
    if ( initialized == false ) {
        import_array();
        initialized = true;
    }
    return initialized;
}
static bool numpy_initialized = _tablerow_initialize_numpy( );


namespace casac {

    // constructor used by from python to construct a tablerow object
    tablerow::tablerow( const casac::table *_table, const std::vector<std::string> &_columnnames, bool _exclude ) :
        itsLog(new casacore::LogIO)
    {
        if ( ! _table ||
             ! _table->itsTable ||
             ! _table->itsTable->isReadable( ) ) throw AipsError( "invalid table passed for parameter one" );
        itsTable = (table*) _table;
        itsProxy = _table->itsTable;
        itsRow.reset( new TableRowProxy( *itsProxy, static_cast<casacore::Vector<casacore::String>>(_columnnames), _exclude ) );
    }

    // constructor used by table class (in table_cmpt.cc) to return a
    // tablerow for fetching one or more rows
    tablerow::tablerow( table *tb, std::shared_ptr<TableHandle> myTable,
                        const std::vector<std::string> &columnnames, bool exclude ) :
        itsLog(new casacore::LogIO), itsProxy(myTable), itsTable(tb)
    {
        if ( ! tb ||
             ! tb->itsTable ||
             ! tb->itsTable->isReadable( ) ) throw AipsError( "invalid table passed for parameter one" );
        itsRow.reset( new TableRowProxy( *itsProxy, static_cast<casacore::Vector<casacore::String>>(columnnames), exclude ) );
    }

    // check to see if tablerow can be modified
    bool tablerow::iswritable( ) {
        *itsLog << LogOrigin(__func__,"");
        try {
            if ( itsRow ) return itsRow->isWritable( ) && itsProxy->isWritable( );
        } catch (AipsError x) {
            *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg( ) << LogIO::POST;
            RETHROW(x);
        }
        reset( );
        throw AipsError( "use of uninitialized table row" );
    }
    // magic function which casacore also supplies
    bool tablerow::_iswritable( ) { return iswritable( ); }

    // fetch one row
    record *tablerow::get( long rownr ) {
        *itsLog << LogOrigin(__func__,"");
        try {
            if ( itsRow ) return fromRecord( itsRow->get( rownr ) );
        } catch (AipsError x) {
            *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg( ) << LogIO::POST;
            RETHROW(x);
        }
        reset( );
        throw AipsError( "use of uninitialized table row" );
    }
    // magic function which casacore also supplies
    record *tablerow::_get( long rownr ) { return get(rownr); }

    // replace one row
    bool tablerow::put( long rownr, const record &value, bool matchingfields) {
        *itsLog << LogOrigin(__func__,"");
        try {
            if ( itsRow ) {
                std::unique_ptr<Record> tab(toRecord( value ));
                itsRow->put( rownr, *tab, matchingfields );
                return true;
            }
        } catch (AipsError x) {
            *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg( ) << LogIO::POST;
            RETHROW(x);
        }
        reset( );
        throw AipsError( "use of uninitialized table row" );
    }
    // magic function casacore also supplies
    bool tablerow::_put( long rownr, const record &value, bool matchingfields) { return put( rownr, value, matchingfields ); }

    // magic function for checking the length (it is not completely clear in which
    // contexts this function is used)
    long tablerow::__len__( ) {
        *itsLog << LogOrigin(__func__,"");
        try {
            if ( itsProxy ) return itsProxy->nrows( );
        } catch (AipsError x) {
            *itsLog << LogIO::SEVERE << "Exception Reported: " << x.getMesg( ) << LogIO::POST;
            RETHROW(x);
        }
        reset( );
        throw AipsError( "use of uninitialized table row" );
    }

    // RAII for PyGILState_Ensure()/Release()
    class GILState_Ensurer {
        PyGILState_STATE state;
        bool inited;
    public:
        void release() { if (inited) { PyGILState_Release(state); inited = false; } }
        GILState_Ensurer() : inited(true), state(PyGILState_Ensure()) {}
        ~GILState_Ensurer() { release(); }
    };

    PyObj* tablerow::__getitem__( PyObj *rownr ) {
        PyObject *obj = (PyObject*) rownr;
        if ( PyNumber_Check(obj) ) {
            // index indicates a single row
            if ( itsProxy && itsRow ) {
                auto pylong = PyNumber_Long(obj);
                auto index = PyLong_AsLong(pylong);
                Py_DECREF(pylong);
                if ( index >= 0 && index < itsProxy->nrows( ) )
                    return toPy( itsRow->get( index ) );
                else
                    throw PyExc_IndexError;
            }
            throw PyExc_IndexError;
        } else if ( PySlice_Check(obj) ) {

            // index indicates a slice
            if ( itsProxy && itsRow ) {
                Py_ssize_t start, stop, step;
                GILState_Ensurer gilState;
                if ( PySlice_Unpack( obj, &start, &stop, &step ) < 0 ) {
                    throw PyExc_IndexError;
                }
                auto slice_length = PySlice_AdjustIndices( itsProxy->nrows( ), &start, &stop, step );
                auto result = PyList_New( slice_length );
                gilState.release();
                for ( ssize_t i=0, row=start; i < slice_length; ++i, row += step ) {
                    if ( row < 0 || row >= itsProxy->nrows( ) ) throw PyExc_IndexError;
                    PyObject *newobj = 0;
                    PyGILState_STATE state;  // Needed for PyGILState_Ensure() and PyGILState_Release()
                    state = PyGILState_Ensure( );
                    try {
                        newobj = toPy( itsRow->get( row ) );
                        PyGILState_Release(state);
                    } catch (...) {
                        PyGILState_Release(state);
                    }
                    if ( PyList_SetItem( result, i, newobj ) < 0 ) {
                        Py_DECREF(result);
                        Py_DECREF(newobj);
                        throw PyExc_ValueError;
                    }
                }
                return result;
            }
            throw PyExc_IndexError;
        } else {
            throw PyExc_IndexError;
        }
        return 0;
    }

    tablerow::~tablerow( ) {
        if ( itsTable ) itsTable->remove_tablerow(this);
        itsProxy.reset( );
        itsTable = 0;
    }

    void tablerow::done( ) {
        if ( itsTable ) itsTable->remove_tablerow(this);
        itsRow.reset( );
        itsProxy.reset( );
        itsTable = 0;
    }

    void tablerow::reset( ) {
        itsRow.reset( );
        itsProxy.reset( );
        itsTable = 0;
    }

}
