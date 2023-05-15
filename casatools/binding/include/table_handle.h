//# table_handle.h
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
//#        Internet email: aips2-request@nrao.edu.
//#        Postal address: AIPS++ Project Office
//#                        National Radio Astronomy Observatory
//#                        520 Edgemont Road
//#                        Charlottesville, VA 22903-2475 USA
//#
#ifndef __CASAC_TABLE_HANDLE_H__
#define __CASAC_TABLE_HANDLE_H__
#include <list>
#include <string>
#include <casacore/casa/aipsxtype.h>
#include <casacore/tables/Tables/TableProxy.h>

using casacore::rownr_t;

namespace casac {

    // The class derived from TableProxy is used to extend the python table binding with
    // additions unique to casatools. To allow extension, this depends upon the casacore
    // PR #1291, and it is used to implement table column iteration for the getcoliter
    // method of the casatools table tool.
    //
    // All of the constructors just forward to the TableProxy base class, and this class
    // does not expand the memory used by TableProxy (so there is no danger of object
    // slicing). It was introduced so that the TableProxy::getValueFromTable(..) function
    // (which is not public) can be reused to implement filling a list with values from
    // the table.
    class TableHandle : public casacore::TableProxy {
      public:
      TableHandle( ) : TableProxy( ) { }
      TableHandle( const TableHandle &tb ) : TableProxy(tb) { }
      TableHandle( const TableProxy &tb ) : TableProxy(tb) { }
      TableHandle( const casacore::Table& table ) : TableProxy(table) { }
      TableHandle( const std::string &tableName,
                   const casacore::Record &lockOptions,
                   int option ) : TableProxy( tableName, lockOptions, option ) { }
      TableHandle( const std::string &tableName,
                   const casacore::Record &lockOptions,
                   const std::string &endianFormat,
                   const std::string &memType,
                   rownr_t nrow,
                   const casacore::Record& tableDesc,
                   const casacore::Record& dmInfo ) 
          : TableProxy( tableName, lockOptions, endianFormat, memType, nrow, tableDesc, dmInfo ) { }
      TableHandle( const std::string &fileName,
                   const std::string &headerName,
                   const std::string &tableName,
                   bool autoHeader,
                   const casacore::IPosition &autoShape,
                   const std::string &separator,
                   const std::string &commentMarker,
                   rownr_t firstLine,
                   rownr_t lastLine,
                   const std::vector<std::string> &columnNames = std::vector<std::string>(),
                   const std::vector<std::string> &dataTypes = std::vector<std::string>() ) 
          : TableProxy( fileName, headerName, tableName, autoHeader, autoShape, separator,
                        commentMarker, firstLine, lastLine,
                        (casacore::Vector<casacore::String>) columnNames,
                        (casacore::Vector<casacore::String>) dataTypes ) { }

        // Read num_rows values from the table column named columnName, starting at start_row
        // with a row increment of incr.
        void fillColumnValues( std::list<casacore::ValueHolder> &dest,
                               const std::string &columnName,
                               rownr_t start_row,
                               rownr_t num_rows,
                               ssize_t incr );
    };

}

#endif
