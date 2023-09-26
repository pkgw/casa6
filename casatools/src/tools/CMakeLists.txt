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

execute_process(COMMAND cat ${CMAKE_CURRENT_SOURCE_DIR}/../../version.txt
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE CASATOOLS_VERSION)
string(REGEX MATCH "([0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+)" CASATOOLS_VERSION ${CASATOOLS_VERSION})

project(
  casatools
  VERSION ${CASATOOLS_VERSION}
  DESCRIPTION "CASATOOLS SWIG wrappers"
  HOMEPAGE_URL https://casa.nrao.edu/
  LANGUAGES CXX)

# Avoid warning messages about old SWIG target name conventions
cmake_policy(SET CMP0078 NEW)
cmake_policy(SET CMP0086 NEW)

include(GNUInstallDirs)
include(CheckCXXCompilerFlag)
include(FindPkgConfig)

# Set a default build type if none was specified
set(default_build_type "RelWithDebInfo")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build." FORCE)
endif()

# Set C++ standard to c++ 14 on M1 Macs.
# Use c++11 on other platforms..

IF(${CMAKE_SYSTEM_PROCESSOR} MATCHES "arm")
  set(CMAKE_CXX_STANDARD 14)
ELSE()
  set(CMAKE_CXX_STANDARD 11)
ENDIF()

set(CMAKE_CXX_STANDARD_REQUIRED ON)

## Set definitions global to the whole project
add_definitions(-DCASACORE_NEEDS_RETHROW)
add_definitions(-DAIPS_DEBUG)

## Check for dependencies

# Python development
find_package (Python3 REQUIRED COMPONENTS Development NumPy)

# SWIG
find_package(SWIG REQUIRED COMPONENTS python)
include(${SWIG_USE_FILE})

# Xerces dependency
pkg_check_modules(XERCES REQUIRED IMPORTED_TARGET xerces-c)

# casacpp dependency
pkg_check_modules(CASACPP REQUIRED IMPORTED_TARGET casacpp)

# some tools depend on casacore directly
pkg_check_modules(CASACORE REQUIRED IMPORTED_TARGET casacore)

# OpenMP dependency
find_package(OpenMP)

# Set version
set(CASATOOLS_MAJOR "${CMAKE_PROJECT_VERSION_MAJOR}")
set(CASATOOLS_MINOR "${CMAKE_PROJECT_VERSION_MINOR}")
set(CASATOOLS_PATCH "${CMAKE_PROJECT_VERSION_PATCH}")
set(CASATOOLS_FEATURE "${CMAKE_PROJECT_VERSION_TWEAK}")
set(CASA_VERSION_STRING "${CMAKE_PROJECT_VERSION}")
configure_file("${CMAKE_SOURCE_DIR}/utils/toolversion.h.in" "${CMAKE_SOURCE_DIR}/utils/toolversion.h" @ONLY)

# Rpath option
option(ENABLE_RPATH "Include rpath in shared libraries and executables" YES)
if(ENABLE_RPATH)
  # Optional but helpful for macosx (install_name field)
  set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}" )
  # Ensure (but this is expected default:)
  set(CMAKE_SKIP_BUILD_RPATH FALSE)
  # Use ORIGIN / rpaths within the build area
  # (python build-install path of wheel would be meaningless)
  set(CMAKE_BUILD_RPATH_USE_ORIGIN TRUE)

  # macOS requires this variable to set the dir part of the lib install_name
  # to rpath but then also manual addition of loader_path paths, as
  # CMAKE_BUILD_RPATH_USE_ORIGIN does not do anything on macOS
  set(CMAKE_MACOSX_RPATH TRUE)
  if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    list(APPEND CMAKE_BUILD_RPATH "@loader_path/lib" "@loader_path")
  endif()

  # Also add auto-found dirs (outside build) to install rpath
  set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  # Ensure (but this is expected default:)
  set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
endif()

# Create a list with all the tools
set(casatools_all_components
  synthesisimstore
  imagepol
  agentflagger
  quanta
  atcafiller
  regionmanager
  synthesismaskhandler
  image
  synthesisnormalizer
  iterbotsink
  sakura
  atmosphere
  sdm                   
  synthesisutils
  linearmosaic
  calanalysis
  sidebandseparator
  logsink
  calibrater
  table
  componentlist
  simulator
  measures
  tablerow
  miriadfiller   
  singledishms
  coordsys
  utils
  spectralline
  vlafiller
  functional
  ms
  msuvbinner
  imager
  synthesisdeconvolver
  msmetadata
  vpmanager
  imagemetadata  
  mstransformer
  synthesisimager
)

# For each tool get the sources and adds them to the preliminary casatools sources
# Also add the tool directory to the include list
foreach(_component IN LISTS casatools_all_components)
  set(PRIVATE_INCLUDES ${PRIVATE_INCLUDES} ${CMAKE_SOURCE_DIR}/${_component})
  file(GLOB_RECURSE casatools_component_sources "${_component}/*.cc")
  set(casatools_library_sources ${casatools_library_sources} ${casatools_component_sources})
endforeach()

# Add the casatools library target to cmake
add_library(casatools SHARED)
# Place lib .so inside build dir for packaging/wheel
set_target_properties (casatools PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/casatools/__casac__/lib")
set_property(TARGET casatools PROPERTY SUFFIX ".${Python3_SOABI}${CMAKE_SHARED_LIBRARY_SUFFIX}")

# Specify which files belong to the casatools target
# Add the binding files to the list of casatool files to compile
file(GLOB_RECURSE casatools_binding_sources ${CMAKE_SOURCE_DIR}/../../binding/source/*.cc)
set(casatools_library_sources ${casatools_library_sources} ${casatools_binding_sources})
# Exclude the *_cmpt.cc files from the list of casatools sources
list(FILTER casatools_library_sources EXCLUDE REGEX "_cmpt.cc")
target_sources(casatools
  PRIVATE 
    ${casatools_library_sources})

# Install the casatools library
install(TARGETS casatools
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Include the generated files by the XML parsing to the includes list of casatools
target_include_directories(casatools PRIVATE ${PRIVATE_INCLUDES})

# Include SWIG bindings to the includes list of casatools
target_include_directories(casatools PRIVATE ${CMAKE_SOURCE_DIR}/../../binding/include)

# CASACPP dependency
target_link_libraries(casatools PUBLIC PkgConfig::CASACPP)

# Other casatools external dependencies
target_link_libraries(casatools PRIVATE Python3::NumPy)
target_link_libraries(casatools PUBLIC PkgConfig::XERCES)
target_link_libraries(casatools PUBLIC readline)
# OpenMP must be set globally, since it adds compiler flags that need to be
# set for certain casacore header or template files, 
# even if the component does not use directly OpenMP
if (OPENMP_FOUND)
  target_link_libraries(casatools PUBLIC OpenMP::OpenMP_CXX)
endif()

# Macro to define each tool
macro(declare_casatool_component name)
  # Define swig file as c++. That configures the swig command
  # to be called n C++ mode and generate C++ bindings
  set_property(SOURCE ${name}/${name}.i PROPERTY CPLUSPLUS ON)
  # Add a SWIG library
  # The SWIG library contains _only_ the *_cmpt.cc file, the rest goes in casatools
  set (UseSWIG_TARGET_NAME_PREFERENCE STANDARD)
  swig_add_library(${name} TYPE SHARED LANGUAGE python SOURCES ${name}/${name}.i ${name}/${name}_cmpt.cc
                   OUTPUT_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/casatools/__casac__" OUTFILE_DIR "${CMAKE_CURRENT_BINARY_DIR}")
  set_property(TARGET ${name} PROPERTY SWIG_INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/../../binding/include)
  set_property(TARGET ${name} PROPERTY SWIG_GENERATED_INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/${name} ${CMAKE_SOURCE_DIR}/../../binding/include)
  set_property(TARGET ${name} PROPERTY SWIG_COMPILE_OPTIONS -threads -arch arm64)

  set_property(TARGET ${name} PROPERTY INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/${name} ${CMAKE_SOURCE_DIR}/../../binding/include)

  # Use .so as library extension rather than CMAKE_SHARED_LIBRARY_SUFFIX, etc. for macos
  # (https://bugs.python.org/issue43898)
  set_property(TARGET ${name} PROPERTY SUFFIX ".${Python3_SOABI}.so")
  set_property(TARGET ${name} PROPERTY LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/casatools/__casac__")
  
  # Add the includes for all tools, since there are some dependencies between them
  target_include_directories(${name} PRIVATE ${PRIVATE_INCLUDES})

  # Add dependencies
  target_link_libraries(${name} PRIVATE PkgConfig::CASACPP)
  target_link_libraries(${name} PRIVATE casatools)
  target_link_libraries(${name} PRIVATE Python3::Module)
  target_link_libraries(${name} PRIVATE Python3::NumPy)
  # OpenMP must be set globally, since it adds compiler flags that need to be
  # set for certain casacore header or template files, 
  # even if the component does not use directly OpenMP
  if (OPENMP_FOUND)
    target_link_libraries(${name} PUBLIC OpenMP::OpenMP_CXX)
  endif()
endmacro()


# Define each tool
foreach(_component IN LISTS casatools_all_components)
  declare_casatool_component(${_component})
endforeach()


# Inter-task dependencies
target_link_libraries(table PRIVATE tablerow)
target_link_libraries(ms PRIVATE table)
target_link_libraries(ms PRIVATE msmetadata)
target_link_libraries(synthesisimager PRIVATE synthesisimstore)
target_link_libraries(synthesisnormalizer PRIVATE synthesisimstore)
target_link_libraries(image PRIVATE componentlist)
target_link_libraries(image PRIVATE coordsys)
target_link_libraries(imagepol PRIVATE image)

# Tasks that depend on casacore
target_link_libraries(table PRIVATE tablerow PkgConfig::CASACORE)
target_link_libraries(ms PRIVATE msmetadata PkgConfig::CASACORE)
target_link_libraries(coordsys PRIVATE PkgConfig::CASACORE)
target_link_libraries(utils PRIVATE PkgConfig::CASACORE)

