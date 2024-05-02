//# table_handle.cc
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
#include <iterator>
#include <casacore/tables/Tables/TableColumn.h>
#include <casacore/casa/Containers/ValueHolder.h>
#include <table_handle.h>
#include <iostream>


namespace casac {

    void TableHandle::fillColumnValues( std::list<casacore::ValueHolder> &dest,
                                        const std::string &columnName,
                                        rownr_t row, rownr_t nrow, ssize_t incr ) {
        // check the number of rows requested based upon the starting row and the row increment
        // (getRowsCheck is provided by the TableProxy class)
        casacore::Int64 nrows = getRowsCheck( columnName, row, nrow, incr, "fillColumnValues" );
        casacore::TableColumn tabcol( table_p, columnName );

        for ( casacore::Int64 i=0; i < nrows; i++ ) {
            if ( tabcol.isDefined(row) ) {
                // Add the result to the destination
                // (getValueFromTable is provided by the TableProxy class)
                dest.push_back( getValueFromTable(columnName, row, 1, 1, false) );
            } else {
                // If a cell value is not defined, add an empty array
                dest.push_back( casacore::ValueHolder( casacore::Array<bool>( ) ) );
            }
            row += incr;
        }
    }

}
