#-------------------------------------------------------------------------------
# SPEX/CMakeLists.txt:  root CMake build rules
#-------------------------------------------------------------------------------

# Copyright (c) 2023-2024, Timothy A. Davis, All Rights Reserved. FIXME
# Just this particular file is under the Apache-2.0 license; each package has
# its own license.
# SPDX-License-Identifier: Apache-2.0

# This file and most packages in SuiteSparse require cmake 3.22 or later.

cmake_minimum_required ( VERSION 3.22 )

project ( "SuiteSparse" )

#-------------------------------------------------------------------------------
# SuiteSparse policies
#-------------------------------------------------------------------------------

set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${PROJECT_SOURCE_DIR}/SuiteSparse_config/cmake_modules )

#-------------------------------------------------------------------------------
# build options
#-------------------------------------------------------------------------------

# lower-case list of all projects that can be built by this root CMake file
set ( SUITESPARSE_ALL_PROJECTS
    "suitesparse_config;amd;colamd;spex" )

# lower-case list of extra projects that can be built by this root CMake file
set ( SUITESPARSE_EXTRA_PROJECTS
    "" )

# lower-case list of known projects that can be built by this root CMake file
set ( SUITESPARSE_KNOWN_PROJECTS "${SUITESPARSE_ALL_PROJECTS};${SUITESPARSE_EXTRA_PROJECTS}" )

set ( SUITESPARSE_ENABLE_PROJECTS "all" CACHE STRING
    "Semicolon-separated list of SuiteSparse projects to be built (${SUITESPARSE_KNOWN_PROJECTS}, or \"all\")" )

# expand "all" early on
if ( "all" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
    set ( SUITESPARSE_ENABLE_PROJECTS "${SUITESPARSE_ENABLE_PROJECTS};${SUITESPARSE_ALL_PROJECTS}" )
    list ( REMOVE_ITEM SUITESPARSE_ENABLE_PROJECTS "all" )
    list ( REMOVE_DUPLICATES  SUITESPARSE_ENABLE_PROJECTS )
endif ( )

# check for unknown projects in list
foreach ( proj ${SUITESPARSE_ENABLE_PROJECTS} )
  if ( NOT "${proj}" IN_LIST SUITESPARSE_KNOWN_PROJECTS )
     message ( FATAL_ERROR "${proj} is not a known project: ${SUITESPARSE_KNOWN_PROJECTS}." )
  endif ( )
endforeach ( )

# options to build with libraries installed on the system instead of building
# dependencies automatically
option ( SUITESPARSE_USE_SYSTEM_AMD "ON: use AMD libraries installed on the build system.  OFF (default): Automatically build AMD as dependency if needed." OFF )
option ( SUITESPARSE_USE_SYSTEM_COLAMD "ON: use COLAMD libraries installed on the build system.  OFF (default): Automatically build COLAMD as dependency if needed." OFF )
option ( SUITESPARSE_USE_SYSTEM_SUITESPARSE_CONFIG "ON: use SuiteSparse_config libraries installed on the build system.  OFF (default): Automatically build SuiteSparse_config as dependency if needed." OFF )

#-------------------------------------------------------------------------------
# global variables
#-------------------------------------------------------------------------------

# Set to indicate that we are building from a root CMake file.
# That will change the directory layout and (imported) target names (namespace!)
# during the build process.
set ( SUITESPARSE_ROOT_CMAKELISTS ON )

#-------------------------------------------------------------------------------
# common SuiteSparse modules
#-------------------------------------------------------------------------------

set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/SuiteSparse_config/cmake_modules )

include ( SuiteSparsePolicy )

#-------------------------------------------------------------------------------
# check/add project dependencies
#-------------------------------------------------------------------------------

if ( SUITESPARSE_USE_SYSTEM_AMD )
    list ( REMOVE_ITEM SUITESPARSE_ENABLE_PROJECTS "amd" )
    find_package ( AMD 3.3.2 REQUIRED )
else ( )
    if ( NOT "amd" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
        message ( STATUS "Adding \"amd\" to the list of built targets." )
        list ( APPEND SUITESPARSE_ENABLE_PROJECTS "amd" )
    endif ( )
endif ( )

if ( SUITESPARSE_USE_SYSTEM_COLAMD )
    list ( REMOVE_ITEM SUITESPARSE_ENABLE_PROJECTS "colamd" )
    find_package ( COLAMD 3.3.3 REQUIRED )
else ( )
    if ( NOT "colamd" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
        message ( STATUS "Adding \"colamd\" to the list of built targets." )
        list ( APPEND SUITESPARSE_ENABLE_PROJECTS "colamd" )
    endif ( )
endif ( )

if ( SUITESPARSE_USE_SYSTEM_SUITESPARSE_CONFIG )
    list ( REMOVE_ITEM SUITESPARSE_ENABLE_PROJECTS "suitesparse_config" )
    find_package ( SuiteSparse_config 7.7.0 REQUIRED )
else ( )
    if ( NOT "suitesparse_config" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
        message ( STATUS "Adding \"suitesparse_config\" to the list of built targets." )
        list ( APPEND SUITESPARSE_ENABLE_PROJECTS "suitesparse_config" )
    endif ( )
endif ( )

if ( CMAKE_VERSION VERSION_LESS 3.24 )
    # work around missing GLOBAL option of find_package in older CMake versions
    # If SuiteSparse is included as a sub-project in other projects, they might
    # need to manually import the OpenMP targets for older CMake versions, too.
    find_package ( OpenMP COMPONENTS C )
endif ( )


#-------------------------------------------------------------------------------
# include selected projects
#-------------------------------------------------------------------------------

if ( "suitesparse_config" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
    add_subdirectory ( "SuiteSparse_config" )
    if ( TARGET SuiteSparseConfig )
        add_library ( SuiteSparse::SuiteSparseConfig ALIAS SuiteSparseConfig )
    else ( )
        add_library ( SuiteSparse::SuiteSparseConfig ALIAS SuiteSparseConfig_static )
    endif ( )
    if ( TARGET SuiteSparseConfig_static )
        add_library ( SuiteSparse::SuiteSparseConfig_static ALIAS SuiteSparseConfig_static )
    endif ( )
endif ( )

if ( "amd" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
    add_subdirectory ( "AMD" )
    if ( TARGET AMD )
        add_library ( SuiteSparse::AMD ALIAS AMD )
    else ( )
        add_library ( SuiteSparse::AMD ALIAS AMD_static )
    endif ( )
    if ( TARGET AMD_static )
        add_library ( SuiteSparse::AMD_static ALIAS AMD_static )
    endif ( )
endif ( )

if ( "colamd" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
    add_subdirectory ( "COLAMD" )
    if ( TARGET COLAMD )
        add_library ( SuiteSparse::COLAMD ALIAS COLAMD )
    else ( )
        add_library ( SuiteSparse::COLAMD ALIAS COLAMD_static )
    endif ( )
    if ( TARGET COLAMD_static )
        add_library ( SuiteSparse::COLAMD_static ALIAS COLAMD_static )
    endif ( )
endif ( )

if ( "spex" IN_LIST SUITESPARSE_ENABLE_PROJECTS )
    add_subdirectory ( "SPEX" )
    if ( TARGET SPEX )
        add_library ( SuiteSparse::SPEX ALIAS SPEX )
    else ( )
        add_library ( SuiteSparse::SPEX ALIAS SPEX_static )
    endif ( )
    if ( TARGET SPEX_static )
        add_library ( SuiteSparse::SPEX_static ALIAS SPEX_static )
    endif ( )
endif ( )

#-------------------------------------------------------------------------------
# report status
#-------------------------------------------------------------------------------

include ( SuiteSparseReport )

#-------------------------------------------------------------------------------
# enable testing facilities
#-------------------------------------------------------------------------------

# Currently, only LAGraph, Mongoose, and CHOLMOD have ctests.

# FIXME: convert more of the existing demos to ctests.

# Most packages have a ./Tcov folder with a full statement coverage test,
# but these are not imported into cmake yet.

# Most packages also have a ./Demo folder, with shorter examples.  These would
# be nice to add as quick ctests.

# CHOLMOD/Tcov takes about 20 minutes to run.  It is also a full coverage
# test of AMD, CAMD, COLAMD, and CCOLAMD, however.  The current CHOLMOD
# ctest is based on a few ./Demo programs.  It's fast but not a full coverate
# test.

# The CSparse/CXSparse Tcov tests are very fast and would be good candidates to
# add.

include ( CTest )

#-------------------------------------------------------------------------------
# rule to remove all files in build directory
#-------------------------------------------------------------------------------

file ( GLOB SUITESPARSE_BUILT_FILES ${PROJECT_BINARY_DIR}/* )
file ( REAL_PATH ${PROJECT_SOURCE_DIR} _real_project_source_dir )
file ( REAL_PATH ${PROJECT_BINARY_DIR} _real_project_binary_dir )
if ( _real_project_source_dir STREQUAL _real_project_binary_dir )
    add_custom_target ( purge
        COMMENT "The target 'purge' is a no-op for in-tree builds.  Consider building out of the source tree." )
else ( )
    add_custom_target ( purge
        COMMAND ${CMAKE_COMMAND} -E echo "Removing files..."
        COMMAND ${CMAKE_COMMAND} -E rm -rf ${SUITESPARSE_BUILT_FILES}
        COMMENT "Purge all files in the build tree" )
endif ( )
