//# SDMSManager.cc: this defines single dish MS transform manager
//#                inheriting MSTransformManager.
//#
//# Copyright (C) 2015
//# National Astronomical Observatory of Japan
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
//# $Id$
#include <algorithm>
#include <array>
#include <iostream>
#include <iterator>
#include <list>

//#include <libsakura/sakura.h>
//#include <libsakura/config.h>

#include <casacore/casa/Logging/LogIO.h>
#include <casacore/casa/Logging/LogOrigin.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Utilities/Sort.h>

#include <casacore/ms/MSSel/MSSelectionTools.h>
#include <msvis/MSVis/VisibilityIterator2.h>
#include <msvis/MSVis/VisSetUtil.h>

#include <casa_sakura/SakuraUtils.h>
#include <singledish/SingleDish/SDMSManager.h>


#define _ORIGIN LogOrigin("SDMSManager", __func__, WHERE)

using namespace casacore;

namespace casa {

SDMSManager::SDMSManager() : doSmoothing_(false) {
}

// SDMSManager &SDMSManager::operator=(SDMSManager const &other)
// {
//   return *this;
// }

SDMSManager::~SDMSManager() {
  LogIO os(_ORIGIN);
}

// -----------------------------------------------------------------------
// Fill output MS with data from an input VisBuffer
// -----------------------------------------------------------------------
void SDMSManager::fillCubeToOutputMs(vi::VisBuffer2 *vb,Cube<Float> const &data_cube,
                                     Cube<Bool> const *flag_cube, Matrix<Float> const *weight_matrix) {
  setupBufferTransformations(vb);

  // make sure the shape of cube matches the number of rows to add.
  AlwaysAssert(data_cube.nplane() == nRowsToAdd_p, AipsError);

  if (bufferMode_p) {
    return;
  }

  // Create RowRef object to fill new rows
  uInt currentRows = outputMs_p->nrow();
  RefRows rowRef(currentRows, currentRows + nRowsToAdd_p / nspws_p - 1);

  // Add new rows to output MS
  outputMs_p->addRow(nRowsToAdd_p, false);

  // Fill new rows
  if (weight_matrix == nullptr) {
    weightSpectrumFlatFilled_p = false;
    weightSpectrumFromSigmaFilled_p = false;
    fillWeightCols(vb, rowRef);
  }
  fillCubeToDataCols(vb, rowRef, data_cube, flag_cube);
  fillIdCols(vb, rowRef);
  if (weight_matrix != nullptr) { // for update_weight=True (CAS-13161)
    outputMsCols_p->weight().putColumnCells(rowRef, *weight_matrix);
  }
}

// ----------------------------------------------------------------------------------------
// Fill main (data) columns which have to be combined together to produce bigger SPWs
// ----------------------------------------------------------------------------------------
void SDMSManager::fillCubeToDataCols(vi::VisBuffer2 *vb, RefRows &rowRef,
                                     Cube<Float> const &data_cube, Cube<Bool> const *flag_cube) {
  ArrayColumn<Bool> *outputFlagCol = NULL;
  for (dataColMap::iterator iter = dataColMap_p.begin(); iter != dataColMap_p.end(); iter++) {
    // Get applicable *_SPECTRUM (copy constructor uses reference semantics)
    // If channel average or combine, otherwise no need to copy
    const Cube<Float> applicableSpectrum = getApplicableSpectrum(vb, iter->first);

    // Apply transformations
    switch (iter->first) {
    case MS::DATA:
      {
        if (mainColumn_p == MS::DATA) {
          outputFlagCol = &(outputMsCols_p->flag());
          setTileShape(rowRef, outputMsCols_p->flag());
        } else {
          outputFlagCol = NULL;
        }
        setTileShape(rowRef, outputMsCols_p->data());
        Cube<Complex> cdata_cube(data_cube.shape());
        convertArray(cdata_cube, data_cube);
        transformCubeOfData(vb, rowRef, cdata_cube, outputMsCols_p->data(),
                            outputFlagCol, applicableSpectrum);
        break;
      }
    case MS::CORRECTED_DATA:
      {
        if (mainColumn_p == MS::CORRECTED_DATA) {
          outputFlagCol = &(outputMsCols_p->flag());
          setTileShape(rowRef, outputMsCols_p->flag());
        } else {
          outputFlagCol = NULL;
        }
        Cube<Complex> cdata_cube(data_cube.shape());
        convertArray(cdata_cube, data_cube);
        if (iter->second == MS::DATA) {
          setTileShape(rowRef, outputMsCols_p->data());
          transformCubeOfData(vb, rowRef, cdata_cube, outputMsCols_p->data(),
                              outputFlagCol, applicableSpectrum);
        } else {
          setTileShape(rowRef, outputMsCols_p->correctedData());
          transformCubeOfData(vb, rowRef, cdata_cube, outputMsCols_p->correctedData(),
                              outputFlagCol, applicableSpectrum);
        }
        break;
      }
    case MS::MODEL_DATA:
      {
        if (mainColumn_p == MS::MODEL_DATA) {
          outputFlagCol = &(outputMsCols_p->flag());
          setTileShape(rowRef, outputMsCols_p->flag());
        } else {
          outputFlagCol = NULL;
        }

        if (iter->second == MS::DATA) {
          setTileShape(rowRef, outputMsCols_p->data());
          transformCubeOfData(vb, rowRef, vb->visCubeModel(), outputMsCols_p->data(),
                              outputFlagCol, applicableSpectrum);
        } else {
          setTileShape(rowRef, outputMsCols_p->modelData());
          transformCubeOfData(vb, rowRef, vb->visCubeModel(), outputMsCols_p->modelData(),
                              outputFlagCol, applicableSpectrum);
        }
        break;
      }
    case MS::FLOAT_DATA:
      {
        if (mainColumn_p == MS::FLOAT_DATA) {
          outputFlagCol = &(outputMsCols_p->flag());
          setTileShape(rowRef, outputMsCols_p->flag());
        } else {
          outputFlagCol = NULL;
        }
        setTileShape(rowRef, outputMsCols_p->floatData());
        transformCubeOfData(vb, rowRef, data_cube, outputMsCols_p->floatData(),
                            outputFlagCol, applicableSpectrum);
        break;
      }
    case MS::LAG_DATA:
      {
        // jagonzal: TODO
        break;
      }
    default:
      {
        // jagonzal: TODO
        break;
      }
    } // end switch

    //KS: THIS PART ASSUMES INROW==OUTROW
    if ((outputFlagCol != NULL) && (flag_cube != nullptr)) {
      setTileShape(rowRef, outputMsCols_p->flag());
      writeCube(*flag_cube, *outputFlagCol, rowRef);
    }
  } // end loop for dataColMap

  // Special case for flag category
  if (inputFlagCategoryAvailable_p) {
    if (spectrumReshape_p) {
      IPosition transformedCubeShape = getShape(); //[nC,nF,nR]
      IPosition inputFlagCategoryShape = vb->flagCategory().shape(); // [nC,nF,nCategories,nR]
      IPosition flagCategoryShape(4,
                                  inputFlagCategoryShape(1),
                                  transformedCubeShape(2),
                                  inputFlagCategoryShape(2),
                                  transformedCubeShape(2));
      Array<Bool> flagCategory(flagCategoryShape, false);

      outputMsCols_p->flagCategory().putColumnCells(rowRef, flagCategory);
    } else {
      outputMsCols_p->flagCategory().putColumnCells(rowRef, vb->flagCategory());
    }
  }

  return;
}

// -----------------------------------------------------------------------
// Method to determine sort columns order
// -----------------------------------------------------------------------
void SDMSManager::setSortColumns(Block<Int> sortColumns,
                                 bool addDefaultSortCols,
                                 Double timebin) {
  size_t const num_in = sortColumns.nelements();
  LogIO os(_ORIGIN);
  os << "Setting user sort columns with " << num_in << " elements" << LogIO::POST;
  if (addDefaultSortCols) {
    Block<Int> defaultCols(6);
    defaultCols[0] = MS::OBSERVATION_ID;
    defaultCols[1] = MS::ARRAY_ID;
    // Assuming auto-correlation here
    defaultCols[2] = MS::ANTENNA1;
    defaultCols[3] = MS::FEED1;
    defaultCols[4] = MS::DATA_DESC_ID;
    defaultCols[5] = MS::TIME;
    os << "Adding default sort columns to user sort column" << LogIO::POST;
    if (num_in > 0) {
      size_t num_elem(num_in);
      sortColumns.resize(num_in + defaultCols.nelements(), false, true);
      for (size_t i = 0; i < defaultCols.nelements(); ++i) {
        bool addcol = true;
        if (getBlockId(sortColumns, defaultCols[i]) > -1) {
          addcol = false; // the columns is in sortColumns
        }
        if (addcol) {
          sortColumns[num_elem] = defaultCols[i];
          ++num_elem;
        }
      }
      // shrink block
      sortColumns.resize(num_elem, true, true);
      userSortCols_ = sortColumns;
    } else {
      userSortCols_ = defaultCols;
    }
  } else { // do not add
    if (num_in > 0) {
      userSortCols_ = sortColumns;
    } else {
      userSortCols_ = Block<Int>();
    }
  }
  os << "Defined user sort columns with " << userSortCols_.nelements() << " elements" << LogIO::POST;
  timeBin_p = timebin;
  os << "Time bin is " << timeBin_p << " sec" << LogIO::POST;
  // regenerate iterator if necessary
  if (visibilityIterator_p != NULL) {
    os << "Regenerating iterator" << LogIO::POST;
    setIterationApproach();
    generateIterator();
  }

  return;
}

int SDMSManager::getBlockId(Block<Int> const &data, Int const value) {
  for (size_t i = 0; i < data.nelements(); ++i) {
    if (data[i] == value) {
      return (int)i;
    }
  }

  return -1;
}

void SDMSManager::setIterationApproach() {
  if (userSortCols_.nelements() == 0) {
    sortColumns_p = Block<Int>();

    return;
  }

  // User column is set.
  logger_p.origin(_ORIGIN);

  using ColumnId = MSMainEnums::PredefinedColumns;
  struct CheckItem {
    ColumnId columnId;
    String columnName;
    String paramType;
    Bool removeIt;
    Bool addIt;
  };
  auto const userSortColExists = [&](ColumnId const &columnId) {
    return getBlockId(userSortCols_, columnId) > -1;
  };

  // copy userSortCols_ to std::list
  std::list<ColumnId> userSortColsList;
  std::transform(
    userSortCols_.begin(), userSortCols_.end(),
    std::back_inserter(userSortColsList),
    [](Int const &i) {return static_cast<ColumnId>(i);}
  );

  // addIt for DATA_DESC_ID is always false because no add operation is intended
  std::array<CheckItem, 4> const checkList {{
    {MS::SCAN_NUMBER, "SCAN_NUMBER", "scan", timespan_p.contains("scan"), !timespan_p.contains("scan")},
    {MS::STATE_ID, "STATE_ID", "state", timespan_p.contains("state"), !timespan_p.contains("state")},
    {MS::FIELD_ID, "FIELD_ID", "field", timespan_p.contains("field"), !timespan_p.contains("field")},
    {MS::DATA_DESC_ID, "DATA_DESC_ID", "spw", combinespws_p, false}
  }};

  // remove columns from userSortColsList if necessary
  for (auto const &item: checkList) {
    if (item.removeIt && userSortColExists(item.columnId)) {
      logger_p << LogIO::NORMAL;
      if (item.paramType.matches("spw")) {
        logger_p << "Combining data from selected spectral windows. ";
      } else {
        logger_p << "Combining data through " << item.paramType << "s for time average. ";
      }
      logger_p << "Removing " << item.columnName << " from user sort list." << LogIO::POST;
      userSortColsList.remove(item.columnId);
    }
  }

  // add columns to userSortColsList if necessary
  if (timeAverage_p) {
    for (auto const &item: checkList) {
      if (item.addIt && !userSortColExists(item.columnId)) {
        logger_p << LogIO::NORMAL
                 << "Splitting data by " << item.paramType << "s for time average. "
                 << "Adding " << item.columnName << " to user sort list." << LogIO::POST;
        userSortColsList.push_back(item.columnId);
      }
    }
  }

  // copy back userSortColsList to sortColumns_p
  sortColumns_p.resize(userSortColsList.size(), true, false);
  std::transform(
    userSortColsList.begin(), userSortColsList.end(), sortColumns_p.begin(),
    [](ColumnId const &i) {return static_cast<Int>(i);}
  );

  ostringstream oss;
  for (size_t i = 0; i < sortColumns_p.nelements(); ++i) {
    oss << sortColumns_p[i] << ", ";
  }
  logger_p << LogIO::DEBUG1 << "Final Sort order: "
           << oss.str() << LogIO::POST;

  return;
}

Record SDMSManager::getSelRec(string const &spw) {
  MeasurementSet myms = MeasurementSet(inpMsName_p);
  return myms.msseltoindex(toCasaString(spw));
}

/*
MeasurementSet SDMSManager::getMS()
{
  MeasurementSet myms = MeasurementSet(inpMsName_p);
  return myms;
}
*/

void SDMSManager::setSmoothing(string const &kernelType, float const &kernelWidth) {
  // kernel type
  VectorKernel::KernelTypes type = VectorKernel::toKernelType(kernelType);

  // Fail if type is neither GAUSSIAN nor BOXCAR since other ones are yet to be supported
  if ((type != VectorKernel::GAUSSIAN) && (type != VectorKernel::BOXCAR)) {
    stringstream oss;
    oss << "Smoothing kernel type \"" << kernelType << "\" is not supported yet.";
    throw AipsError(oss.str());
  }

  // Fail if kernel width is zero or negative
  if (kernelWidth <= 0.0) {
    throw AipsError("Zero or negative kernel width is not allowed.");
  }

  doSmoothing_ = true;
  kernelType_ = type;
  kernelWidth_ = kernelWidth;
}

void SDMSManager::unsetSmoothing() {
  doSmoothing_ = false;
}

void SDMSManager::initializeSmoothing() {
  LogIO os(_ORIGIN);
  if (!doSmoothing_) {
    return;
  }

  Vector<Int> numChanList = inspectNumChan();
  for (size_t i = 0; i < numChanList.nelements(); ++i) {
    Int numChan = numChanList[i];
    Vector<Float> theKernel = VectorKernel::make(kernelType_, kernelWidth_, numChan, true, false);
    //shift 1 channel for boxcar kernel---(for CAS-7442, 2015/11/18 WK)
    if (kernelType_ == VectorKernel::BOXCAR) {
      for (size_t j = theKernel.nelements()-1; j >= 1; --j) {
        theKernel[j] = theKernel[j-1];
      }
      theKernel[0] = 0.0;
    }
    //------
    convolverPool_[numChan] = Convolver<Float>();
    convolverPool_[numChan].setPsf(theKernel, IPosition(1, numChan));
  }

  // set smoothing kernel for weight
  smoothBin_p = static_cast<uInt>(kernelWidth_ + 0.5f);
  smoothBin_p += (smoothBin_p % 2 == 0) ? 1 : 0; // to make smoothBin_p odd
  uInt numChanMinimum = min(numChanList);
  smoothBin_p = min(smoothBin_p, numChanMinimum - ((numChanMinimum % 2 == 0) ? 1 : 0)); // smoothBin_p < numChanMinimum
  uInt halfWidth = smoothBin_p / 2;
  Sort sort;
  Vector<Float> kernelForMinimumNumChan = convolverPool_[numChanMinimum].getPsf();
  sort.sortKey(kernelForMinimumNumChan.data(), TpFloat, 0, Sort::Descending);
  Vector<uInt> indexArray;
  //uInt indexArrayLength = sort.sort(indexArray, numChanMinimum);
  sort.sort(indexArray, numChanMinimum);
  uInt startChan = indexArray[0] - halfWidth;
  uInt endChan = startChan + smoothBin_p;
  smoothCoeff_p.resize(smoothBin_p, false);
  for (uInt i = startChan, j = 0; i < endChan; ++i, ++j) {
    smoothCoeff_p[j] = kernelForMinimumNumChan[i];
  }
  // normalize smoothCoeff_p
  smoothCoeff_p /= sum(smoothCoeff_p);
  os << LogIO::DEBUGGING << "smoothBin_p = " << smoothBin_p << LogIO::POST;
  os << LogIO::DEBUGGING << "smoothCoeff_p = " << smoothCoeff_p << LogIO::POST;
}

Vector<Int> SDMSManager::inspectNumChan() {
  LogIO os(_ORIGIN);

  if (selectedInputMs_p == NULL) {
    throw AipsError("Input MS is not opened yet.");
  }

  ScalarColumn<Int> col(*selectedInputMs_p, "DATA_DESC_ID");
  Vector<Int> ddIdList = col.getColumn();
  uInt numDDId = GenSort<Int>::sort(ddIdList, Sort::Ascending,
                                    Sort::QuickSort | Sort::NoDuplicates);
  col.attach(selectedInputMs_p->dataDescription(), "SPECTRAL_WINDOW_ID");
  Vector<Int> spwIdList(numDDId);
  for (uInt i = 0; i < numDDId; ++i) {
    spwIdList[i] = col(ddIdList[i]);
  }
  Vector<Int> numChanList(numDDId);
  col.attach(selectedInputMs_p->spectralWindow(), "NUM_CHAN");
  os << LogIO::DEBUGGING << "spwIdList = " << spwIdList << LogIO::POST;
  for (size_t i = 0; i < spwIdList.nelements(); ++i) {
    Int spwId = spwIdList[i];
    Int numChan = col(spwId);
    numChanList[i] = numChan;
    os << LogIO::DEBUGGING << "examine spw " << spwId << ": nchan = " << numChan << LogIO::POST;
    if (numChan == 1) {
      stringstream ss;
      ss << "smooth: Failed due to wrong spw " << i;
      throw AipsError(ss.str());
    }
  }

  uInt numNumChan = GenSort<Int>::sort(numChanList, Sort::Ascending,
                                       Sort::QuickSort | Sort::NoDuplicates);
  return numChanList(Slice(0, numNumChan));
}

}  // End of casa namespace.
