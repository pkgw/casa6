# -*- mode: cmake -*-
# SPDX-License-Identifier: GPL-2.0-or-later
# Copyright (C) 2022 CASA project
#
# This library is free software; you can redistribute it and/or modify it
# under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
# License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Correspondence concerning CASA should be addressed as follows:
#        Internet email: casa-feedback@nrao.edu.

## Set general cmake configuration

# Minimum cmake required version
cmake_minimum_required(VERSION 3.12)

project(
  casacpp
  VERSION 6.5.2.0
  DESCRIPTION "CASA C++ libraries"
  HOMEPAGE_URL https://casa.nrao.edu/
  LANGUAGES CXX Fortran C)

# Add directory with additional macros
set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

include(DeclareCasacppComponent)
include(GNUInstallDirs)
include(CheckCXXCompilerFlag)
include(FindPkgConfig)

# Set a default build type if none was specified
set(default_build_type RelWithDebInfo)

# Set C++ standard to c++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## Check for dependencies

# Find MPI compiler
find_package(MPI)

# casacore dependency
pkg_check_modules(CASACORE REQUIRED IMPORTED_TARGET casacore)

# Find BISON, FLEX
find_package(BISON REQUIRED)
find_package(FLEX REQUIRED)

# protobuf/grpc dependency
find_package(Protobuf REQUIRED)
find_program(_GRPC_CPP_PLUGIN_EXECUTABLE grpc_cpp_plugin)
# Try to detect GRPC using native cmake support (for distributions like Fedora or source compilation)
find_package(gRPC QUIET)
# If that fails, try detecting GRPC with pkgconfig (for distributions like Ubuntu, Rocky Linux)
if(NOT DEFINED gRPC_FOUND OR gRPC_FOUND EQUAL 0)
  pkg_check_modules(GRPC REQUIRED IMPORTED_TARGET grpc++)
  set(gRPC_USES_FIND_PACKAGE OFF)
else()
  set(gRPC_USES_FIND_PACKAGE ON)
endif()
add_definitions(-DUSE_GRPC)

# GSL dependency
find_package(GSL REQUIRED)

# LIBXML dependecy
pkg_check_modules(LibXML++ REQUIRED IMPORTED_TARGET libxml++-2.6)

# LIBXSLT dependecy
pkg_check_modules(XSLT REQUIRED IMPORTED_TARGET libxslt)

# CFITSIO dependecy
pkg_check_modules(CFITSIO REQUIRED IMPORTED_TARGET cfitsio)

# FFTW dependency
pkg_search_module(FFTW REQUIRED IMPORTED_TARGET fftw3)
pkg_search_module(FFTWF REQUIRED IMPORTED_TARGET fftw3f)
find_library(FFTW_DOUBLE_THREADS_LIB
             NAMES "fftw3_threads"
             PATHS ${PKG_FFTW_LIBRARY_DIRSOOO} REQUIRED)
find_library(FFTW_FLOAT_THREADS_LIB
             NAMES "fftw3f_threads"
             PATHS ${PKG_FFTWF_LIBRARY_DIRS} REQUIRED)

# sqllite3 dependency
find_package(SQLite3 REQUIRED)

# OpenMP dependency
find_package(OpenMP)

# SAKURA dependency
pkg_check_modules(SAKURA REQUIRED IMPORTED_TARGET libsakura)

# Create a list with all the components
set(casacpp_all_components
  alma
  asdmstman
  atmosphere
  atnf
  calanalysis
  casadbus
  casatools
  casa_sakura
  components
  imageanalysis
  flagging
  graphics
  miriad
  mstransform
  msvis
  nrao
  parallel
  protobuf_generated
  singledish
  singledishfiller
  spectrallines
  stdcasa
  synthesis
)

# Add a subdirectory per casacpp component
foreach(_component IN LISTS casacpp_all_components)
  add_subdirectory(${CMAKE_SOURCE_DIR}/${_component})
  set(PRIVATE_LIBS "${PRIVATE_LIBS} -lcasacpp_${_component}")
endforeach()
