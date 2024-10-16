//# table_getcoliter.cc
//# Copyright (C) 2023
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
#include <cmath>
#include <memory>
#include <iostream>
#include <algorithm>
#include <Python.h>
#include <table_cmpt.h>
#include <table_handle.h>
#include <table_convert.h>
#include <casacore/casa/Containers/ValueHolder.h>

// ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
// initialize numpy arrays; this initialization must be done once
// for each c++ tranlation unit.
// ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
static bool _table_getcoliter_initialize_numpy( ) {
    static bool initialized = false;
    if ( initialized == false ) {
        import_array();
        initialized = true;
    }
    return initialized;
}
static bool numpy_initialized = _table_getcoliter_initialize_numpy( );

namespace casac {

    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    // python object storage definition for table column iteration
    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    typedef struct {
        PyObject_HEAD
        bool to_record;                      // should dictionaries be used to return values
        ssize_t incr;                        // row increment when loading from the table
        rownr_t start_row;                   // row where iteration should start
        rownr_t total_to_return;             // the number of rows which shold be returned
        rownr_t total_sent;                  // how many rows have already been sent
        bool iteration_overrun;              // repeated calls to get iteration elements adds error message
        std::shared_ptr<TableHandle>  table; // TableHandle for access to the casacore Table

        //--- lists resize themselves (apparently) so lists cannot be store directly ---
        // cache is a list of cache lists, one per column to be returned
        std::shared_ptr<std::list<std::list<casacore::ValueHolder>>> cache;
        // column_names is a list of the column names
        std::shared_ptr<std::list<std::string>> column_names;
    } getcoliter_Iter;


    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    // retrieving an iterator (from our iterator) just returns a reference to itself
    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    PyObject* getcoliter_Iter_iter(PyObject *self) {
        Py_INCREF(self);
        return self;
    }


    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    // construct a python tuple from the values of one or more column
    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    static PyObject *generate_tuple( std::list<std::list<casacore::ValueHolder>> &values ) {
        if ( values.size( ) == 0 ) {
            // no values are available, this is probably an error
            Py_RETURN_NONE;

        } else if ( values.size( ) == 1 ) {
            // because there is only one column, return the single value not a
            // single element tuple because the is the default behavior with
            // tuples within python
            std::list<casacore::ValueHolder> &lst = values.front( );
            if ( lst.size( ) > 0 ) {
                casacore::ValueHolder vh = lst.front( );
                lst.pop_front( );
                PyObject * result = NULL;
                {
                    PyGILState_STATE state = PyGILState_Ensure( );
                    try {
                        result = toPy(vh);
                        PyGILState_Release(state);
                    } catch (...) {
                        PyGILState_Release(state);
                    }
                }
                if ( result ) return result;
                else Py_RETURN_NONE;
            } else {
                Py_RETURN_NONE;
            }

        } else {
            // with a number of columns return tuple with a set of values
            PyGILState_STATE state = PyGILState_Ensure( );
            PyObject *result = PyTuple_New(values.size( ));
            if ( result == nullptr ) {
                result = PyErr_NoMemory();
                PyGILState_Release(state);
                return result;
            }

            Py_ssize_t index=0;
            for ( auto ptr=values.begin( ); ptr != values.end( ); ++ptr ) {
                casacore::ValueHolder vh = ptr->front( );
                ptr->pop_front( );
                PyObject *pyval = NULL;
                try {
                    pyval = toPy(vh);
                } catch (...) {
                    Py_DECREF(result);
                    PyErr_SetString(PyExc_RuntimeError, "failed to set tuple values" );
                    PyGILState_Release(state);
                    return NULL;
                }
                PyTuple_SetItem( result, index++, pyval );
            }
            PyGILState_Release(state);
            return result;
        }
    }

    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    // construct a python tuple from the values of one or more column
    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    static PyObject *generate_record( std::list<std::string> &keys, std::list<std::list<casacore::ValueHolder>> &values ) {
        if ( values.size( ) == 0 || keys.size( ) == 0 ) {
            // no values are available, this is probably an error
            Py_RETURN_NONE;

        } else {
            // convert available values and names to dictionary elements
            PyGILState_STATE state = PyGILState_Ensure( );
            PyObject *result = PyDict_New( );
            if ( result == nullptr ) {
                result = PyErr_NoMemory();
                PyGILState_Release(state);
                return result;
            }

            auto vptr = values.begin( );
            auto kptr = keys.begin( );
            for ( ; vptr != values.end( ) && kptr != keys.end( ); ++vptr, ++kptr ) {
                casacore::ValueHolder vh = vptr->front( );
                vptr->pop_front( );
                PyObject *value = NULL;
                try {
                    value = toPy(vh);
                } catch (...) {
                    Py_DECREF(result);
                    PyErr_SetString(PyExc_RuntimeError, "failed to set record values" );
                    PyGILState_Release(state);
                    return NULL;
                }
                PyDict_SetItemString( result, kptr->c_str( ), value );
                Py_DECREF( value );
            }
            PyGILState_Release(state);
            return result;
        }
    }

    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    // colum iterator next( ) function which returns the next element from the column
    // (or columns). the return value can be either a tuple or a dictionary
    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    PyObject* getcoliter_Iter_iternext(PyObject *self) {
        constexpr rownr_t ROW_CACHE_SIZE = 512;
        getcoliter_Iter *p = (getcoliter_Iter *)self;
        if ( p->total_sent >= p->total_to_return ) {
            /* Raising of standard StopIteration exception with empty value. */
            if ( p->iteration_overrun )  {
                PyGILState_STATE state = PyGILState_Ensure( );
                PyErr_SetString( PyExc_StopIteration, "attempted iteration beyond the end of iterator" );
                PyGILState_Release(state);
            } else PyErr_SetNone( PyExc_StopIteration );
            p->iteration_overrun = true;
            for ( auto ptr=p->cache->begin( ); ptr != p->cache->end( ); ++ptr )
                ptr->clear( );
            return NULL;
        }

        // check to ensure that a value is available for all columns
        if ( all_of( p->cache->begin( ), p->cache->end( ),
                     [](const std::list<casacore::ValueHolder> &l) { return l.size( ) > 0; } ) ) {
            auto result = p->to_record ? generate_record( *p->column_names, *p->cache ) : generate_tuple( *p->cache );
            (p->total_sent)++;
            return result;
        }

        // at this point at least one of the lists in the cache is zero and none of them
        // are expected to be zero, but if they are, clear them to start fresh
        for ( auto ptr=p->cache->begin( ); ptr != p->cache->end( ); ++ptr )
            ptr->clear( );

        // refill the cache of column values
        auto cache = p->cache->begin( );
        auto name = p->column_names->begin( );
        auto cache_size = std::min( (rownr_t) ceil(ROW_CACHE_SIZE / p->cache->size( )),
                                    p->total_to_return - p->total_sent );
        for ( ; cache != p->cache->end( ) && name != p->column_names->end( ); ++cache, ++name ) {
            try {
                p->table->fillColumnValues( *cache,
                                            *name,
                                            p->start_row + p->total_sent,
                                            cache_size,
                                            p->incr );
            } catch ( const casacore::AipsError &ae ) {
                PyErr_SetString(PyExc_RuntimeError, ae.what( ));
                return NULL;
            }
        }

        // verify that we were able to fill all column cache values
        if ( ! all_of( p->cache->begin( ), p->cache->end( ),
                       [=](const std::list<casacore::ValueHolder> &l) { return l.size( ) == cache_size; } ) ) {
            p->total_to_return = p->total_sent;
            PyGILState_STATE state = PyGILState_Ensure( );
            PyErr_SetString(PyExc_RuntimeError, "loading values from table failed" );
            PyGILState_Release(state);
            return NULL;
        }

        // return the result as a tuple or a record as indicated
        auto result = p->to_record ? generate_record( *p->column_names, *p->cache ) : generate_tuple( *p->cache );
        if ( result ) (p->total_sent)++;
        return result;
    }

    // forward declare rich comparison function because it uses the address of
    // getcoliter_IterType and getcoliter_IterType includes a pointer to
    // getcoliter_richcmp
    static PyObject *getcoliter_richcmp(PyObject *obj1, PyObject *obj2, int op);

    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    // clean up the memory and references when the iteration object is deleted
    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    void getcoliter_IterDealloc( PyObject *self ) {
        PyTypeObject *tp = Py_TYPE(self);
        getcoliter_Iter *myObj = (getcoliter_Iter*) self;
        myObj->table.reset( );
        myObj->cache->clear( );
        myObj->column_names->clear( );
        PyObject_Del( self );
    }

    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    // python type spectification for table column iteration
    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    static PyTypeObject getcoliter_IterType = {
        PyVarObject_HEAD_INIT(&PyType_Type, 0)
        "table._getcoliter_iter",                   /* tp_name */
        sizeof(getcoliter_Iter),                    /* tp_basicsize */
        0,                                          /* tp_itemsize */
        getcoliter_IterDealloc,                     /* tp_dealloc */
        0,                                          /* tp_print */
        0,                                          /* tp_getattr */
        0,                                          /* tp_setattr */
        0,                                          /* tp_compare (python2) or tp_reserved (python3) */
        0,                                          /* tp_repr */
        0,                                          /* tp_as_number */
        0,                                          /* tp_as_sequence */
        0,                                          /* tp_as_mapping */
        0,                                          /* tp_hash  */
        0,                                          /* tp_call */
        0,                                          /* tp_str */
        PyObject_GenericGetAttr,                    /* tp_getattro */
        0,                                          /* tp_setattro */
        0,                                          /* tp_as_buffer */
        Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags: Py_TPFLAGS_HAVE_ITER tells python to
                                                       use tp_iter and tp_iternext fields. */
        "Internal myiter iterator object.",         /* tp_doc */
        0,                                          /* tp_traverse */
        0,                                          /* tp_clear */
        getcoliter_richcmp,                         /* tp_richcompare */
        0,                                          /* tp_weaklistoffset */
        getcoliter_Iter_iter,                       /* tp_iter: __iter__() method */
        getcoliter_Iter_iternext                    /* tp_iternext: next() method */
                                                    /* tp_methods */
                                                    /* tp_members */
                                                    /* tp_getset */
                                                    /* tp_base */
                                                    /* tp_dict */
                                                    /* tp_descr_get */
                                                    /* tp_descr_set */
                                                    /* tp_dictoffset */
                                                    /* tp_init */
                                                    /* tp_alloc */
                                                    /* tp_new */
                                                    /* tp_free */
    };

    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    // this is called for things like:
    //
    //     tb.getcol( 'TIME', 0, 7, 100 ) == tb.getcoliter( 'TIME', 0, 7, 100 )
    //
    // where something like a numpy array is compared to this iterator.
    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    static PyObject *getcoliter_richcmp(PyObject *obj1, PyObject *obj2, int op) {
        // expect obj1 to be a table column iterator
        if ( Py_TYPE(obj1) != &getcoliter_IterType )
            Py_RETURN_NOTIMPLEMENTED;

        // get the next iteration value
        PyObject *next_object = getcoliter_Iter_iternext(obj1);

        // if we're at the end of the iterator, return NULL
        // (and assume StopIteration has been set)...
        if ( next_object == NULL ) { return NULL; }

        // compare the next value against obj2
        PyObject *result = PyObject_RichCompare( next_object, obj2, op );

        // free the retrieved object
        Py_DECREF(next_object);

        return result;
    }

    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    // table object member function which returns the iterator
    // ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----  ----
    PyObj* table::getcoliter( const vector<string>& _columnnames, long _startrow, long _nrow, long _rowincr, bool _torecord ) {

        // if the user has not opened the table, this function can still be called
        // but the itsTable TableProxy member will be NULL
        if ( itsTable == nullptr ) {
            PyGILState_STATE state = PyGILState_Ensure( );
            PyErr_SetString(PyExc_RuntimeError, "no opened table available" );
            PyGILState_Release(state);
            return NULL;
        }

        auto column_names = itsTable->columnNames( );

        // the GIL must be locked to allocate the new object
        getcoliter_Iter *p = 0;
        {
            PyGILState_STATE state = PyGILState_Ensure( );
            try {
                // this must be within the GIL lock otherwise a SEGV occurs
                for ( auto name : _columnnames ) {
                    if ( std::find( begin(column_names), end(column_names), casacore::String(name) ) == std::end(column_names) ) {
                        PyErr_Format( PyExc_RuntimeError, "column \"%s\" does not exist", name.c_str( ) );
                        PyGILState_Release(state);
                        return NULL;
                    }
                }
                // create an object
                p = PyObject_New(getcoliter_Iter, &getcoliter_IterType);
                PyGILState_Release(state);
            } catch (...) {
                PyGILState_Release(state);
            }
        }
        if ( ! p ) {
            PyGILState_STATE state = PyGILState_Ensure( );
            PyErr_SetString(PyExc_RuntimeError, "could not create iterator");
            PyGILState_Release(state);
            return NULL;
        }

        // initialize the iterator object using the type specification
        if ( ! PyObject_Init((PyObject *)p, &getcoliter_IterType)) {
            Py_DECREF(p);
            PyGILState_STATE state = PyGILState_Ensure( );
            PyErr_SetString(PyExc_RuntimeError, "could not create iterator");
            PyGILState_Release(state);
            return NULL;
        }

        // because the struct is malloc'ed the objects it contains
        // must be initialized by extraordinary means...
        // uninitialized shared pointers are used to initialize
        // the pointers in the object's struct
        std::shared_ptr<TableHandle> th_init;
        shared_ptr<std::list<std::list<casacore::ValueHolder>>> cache_init;
        shared_ptr<std::list<std::string>> coln_init;
        memcpy( &p->table, &th_init, sizeof(p->table) );
        memcpy( &p->cache, &cache_init, sizeof(p->cache) );
        memcpy( &p->column_names, &coln_init, sizeof(p->column_names) );

        // initialize the column name and cache list, then fill
        // them with names and empty cache lists
        p->cache.reset(new std::list<std::list<casacore::ValueHolder>>( ));
        p->column_names.reset(new std::list<std::string>( ));
        for ( auto cur=_columnnames.begin( ); cur != _columnnames.end( ); ++cur ) {
            p->column_names->push_back(*cur);
            p->cache->push_back(std::list<casacore::ValueHolder>( ));
        }

        // initialize iteration state
        p->start_row = _startrow < 0 ? 0 : _startrow;
        p->total_to_return = _nrow < 0 ? itsTable->nrows( ) : _nrow;
        p->total_sent = 0;
        p->iteration_overrun = false;
        p->incr = _rowincr < 0 ? 1 : _rowincr;
        p->to_record = _torecord;

        // set TableHandle (derived from TableProxy) reference
        p->table = itsTable;
        return (PyObject*) p;
    }
}
