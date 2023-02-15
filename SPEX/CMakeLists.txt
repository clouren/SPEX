#-------------------------------------------------------------------------------
# SuiteSparse/SPEX/CMakeLists.txt:  cmake for SPEX
#-------------------------------------------------------------------------------

# Copyright (c) 2022-2023, Timothy A. Davis.
# All Rights Reserved.
# SPDX-License-Identifier: BSD-3-clause

#-------------------------------------------------------------------------------
# get the version
#-------------------------------------------------------------------------------

cmake_minimum_required ( VERSION 3.19 )

set ( SPEX_DATE "Dec 27, 2022" )
set ( SPEX_VERSION_MAJOR 3 )
set ( SPEX_VERSION_MINOR 0 )
set ( SPEX_VERSION_SUB   0 )

message ( STATUS "Building SPEX version: v"
    ${SPEX_VERSION_MAJOR}.
    ${SPEX_VERSION_MINOR}.
    ${SPEX_VERSION_SUB} " (" ${SPEX_DATE} ")" )

#-------------------------------------------------------------------------------
# SuiteSparse policies
#-------------------------------------------------------------------------------

set ( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake_modules
    ${CMAKE_SOURCE_DIR}/../COLAMD/cmake_modules
    ${CMAKE_SOURCE_DIR}/../AMD/cmake_modules
    ${CMAKE_SOURCE_DIR}/../SuiteSparse_config/cmake_modules )

include ( SuiteSparsePolicy )

#-------------------------------------------------------------------------------
# define the project
#-------------------------------------------------------------------------------

project ( spex
    VERSION "${SPEX_VERSION_MAJOR}.${SPEX_VERSION_MINOR}.${SPEX_VERSION_SUB}"
    LANGUAGES C )

#-------------------------------------------------------------------------------
# find library dependencies
#-------------------------------------------------------------------------------

find_package ( SuiteSparse_config 6.0.2 REQUIRED )
find_package ( COLAMD 3.0.2 REQUIRED )
find_package ( AMD 3.0.2 REQUIRED )
find_package ( GMP 6.1.2 REQUIRED )     # from SPEX/cmake_modules
find_package ( MPFR 4.0.2 REQUIRED )    # from SPEX/cmake_modules

option ( NOPENMP "ON (default): do not use OpenMP.  OFF: use OpenMP" on )
if ( NOPENMP )
    set ( OPENMP_FOUND false )
else ( )
    find_package ( OpenMP )
endif ( )

include ( SuiteSparse__thread )

#-------------------------------------------------------------------------------
# configure files
#-------------------------------------------------------------------------------

configure_file ( "Config/SPEX.h.in"
                 "${PROJECT_SOURCE_DIR}/Include/SPEX.h")
configure_file ( "Config/SPEX_version.tex.in"
                 "${PROJECT_SOURCE_DIR}/Doc/SPEX_version.tex")
configure_file ( "Config/spex_deps.m.in"
                 "${PROJECT_SOURCE_DIR}/MATLAB/spex_deps.m")

#-------------------------------------------------------------------------------
# include directories
#-------------------------------------------------------------------------------

include_directories ( Include
    SPEX_Backslash/Source SPEX_Cholesky/Source SPEX_LU/Source
    SPEX_Update/Source SPEX_Utilities/Source
    ${SUITESPARSE_CONFIG_INCLUDE_DIR}
    ${GMP_INCLUDE_DIR} ${MPFR_INCLUDE_DIR}
    ${AMD_INCLUDE_DIR} ${COLAMD_INCLUDE_DIR} )

# TODO:
#   SPEX_QR/Source

#-------------------------------------------------------------------------------
# compiler flags
#-------------------------------------------------------------------------------

# gcc -O3 complains about SPEX_mpq_init

# if ( "${CMAKE_C_COMPILER_ID}" STREQUAL "GNU" )
#     set ( CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O2" )
# endif ( )

#-------------------------------------------------------------------------------
# dynamic spex library properties
#-------------------------------------------------------------------------------

file ( GLOB SPEX_SOURCES "SPEX*/Source/*.c" )

add_library ( spex SHARED ${SPEX_SOURCES} )

set_target_properties ( spex PROPERTIES
    VERSION ${SPEX_VERSION_MAJOR}.${SPEX_VERSION_MINOR}.${SPEX_VERSION_SUB}
    C_STANDARD_REQUIRED 11
    SOVERSION ${SPEX_VERSION_MAJOR}
    PUBLIC_HEADER "Include/SPEX.h" )

#-------------------------------------------------------------------------------
# static spex library properties
#-------------------------------------------------------------------------------

if ( NOT NSTATIC )
add_library ( spex_static STATIC ${SPEX_SOURCES} )

set_target_properties ( spex_static PROPERTIES
    VERSION ${SPEX_VERSION_MAJOR}.${SPEX_VERSION_MINOR}.${SPEX_VERSION_SUB}
    C_STANDARD_REQUIRED 11
    OUTPUT_NAME spex
    SOVERSION ${SPEX_VERSION_MAJOR} )
endif ( )

#-------------------------------------------------------------------------------
# add the library dependencies
#-------------------------------------------------------------------------------

# suitesparseconfig:
target_link_libraries ( spex PUBLIC ${SUITESPARSE_CONFIG_LIBRARIES} )
if ( NOT NSTATIC )
    target_link_libraries ( spex_static PUBLIC ${SUITESPARSE_CONFIG_LIBRARIES} )
endif ( )

# OpenMP:
if ( OPENMP_FOUND )
    message ( STATUS "OpenMP C libraries:      ${OpenMP_C_LIBRARIES} ")
    message ( STATUS "OpenMP C include:        ${OpenMP_C_INCLUDE_DIRS} ")
    message ( STATUS "OpenMP C flags:          ${OpenMP_C_FLAGS} ")
    target_link_libraries ( spex PUBLIC ${OpenMP_C_LIBRARIES} )
    if ( NOT NSTATIC )
        target_link_libraries ( spex_static PUBLIC ${OpenMP_C_LIBRARIES} )
    endif ( )
    set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} " )
    include_directories ( ${OpenMP_C_INCLUDE_DIRS} )
endif ( )

# AMD:
target_link_libraries ( spex PUBLIC ${AMD_LIBRARIES} )
if ( NOT NSTATIC )
    target_link_libraries ( spex_static PUBLIC ${AMD_STATIC} )
endif ( )

# COLAMD:
target_link_libraries ( spex PUBLIC ${COLAMD_LIBRARIES} )
if ( NOT NSTATIC )
    target_link_libraries ( spex_static PUBLIC ${COLAMD_STATIC} )
endif ( )

# MPFR:
target_link_libraries ( spex PUBLIC ${MPFR_LIBRARIES} )
if ( NOT NSTATIC )
    target_link_libraries ( spex_static PUBLIC ${MPFR_STATIC} )
endif ( )

# GMP:
# must occur after MPFR
target_link_libraries ( spex PUBLIC ${GMP_LIBRARIES} )
if ( NOT NSTATIC )
    target_link_libraries ( spex_static PUBLIC ${GMP_STATIC} )
endif ( )

# libm:
if ( NOT WIN32 )
    target_link_libraries ( spex PUBLIC m )
    if ( NOT NSTATIC )
        target_link_libraries ( spex_static PUBLIC m )
    endif ( )
endif ( )

#-------------------------------------------------------------------------------
# SPEX installation location
#-------------------------------------------------------------------------------

install ( TARGETS spex
    LIBRARY DESTINATION ${SUITESPARSE_LIBDIR}
    ARCHIVE DESTINATION ${SUITESPARSE_LIBDIR}
    RUNTIME DESTINATION ${SUITESPARSE_BINDIR}
    PUBLIC_HEADER DESTINATION ${SUITESPARSE_INCLUDEDIR} )
install ( FILES ${CMAKE_SOURCE_DIR}/cmake_modules/FindSPEX.cmake
    DESTINATION ${SUITESPARSE_LIBDIR}/cmake/SuiteSparse 
    COMPONENT Development )
if ( NOT NSTATIC )
    install ( TARGETS spex_static
        ARCHIVE DESTINATION ${SUITESPARSE_LIBDIR} )
endif ( )

#-------------------------------------------------------------------------------
# python interface
#-------------------------------------------------------------------------------

file ( GLOB SPEX_PYTHON_SOURCES "Python/Source/*.c" )

add_library ( spexpython SHARED ${SPEX_PYTHON_SOURCES} )

set_target_properties ( spexpython PROPERTIES
    VERSION ${SPEX_VERSION_MAJOR}.${SPEX_VERSION_MINOR}.${SPEX_VERSION_SUB}
    C_STANDARD_REQUIRED 11
    SOVERSION ${SPEX_VERSION_MAJOR}
    PUBLIC_HEADER "Python/Source/SPEX_connect.h" )

target_link_libraries ( spexpython PUBLIC spex )

install ( TARGETS spexpython
    LIBRARY DESTINATION ${SUITESPARSE_LIBDIR}
    ARCHIVE DESTINATION ${SUITESPARSE_LIBDIR}
    RUNTIME DESTINATION ${SUITESPARSE_BINDIR}
    PUBLIC_HEADER DESTINATION ${SUITESPARSE_INCLUDEDIR} )

#-------------------------------------------------------------------------------
# Demo library and programs
#-------------------------------------------------------------------------------

option ( DEMO "ON: Build the demo programs.  OFF (default): do not build the demo programs." off )
if ( DEMO )

    #---------------------------------------------------------------------------
    # demo library
    #---------------------------------------------------------------------------

    message ( STATUS "Also compiling the demos in SPEX/Demo" )

    #---------------------------------------------------------------------------
    # Demo programs
    #---------------------------------------------------------------------------

    include_directories ( Demo )
    file ( GLOB SPEX_DEMO_SOURCES "Demo/Utilities/*.c" )

    add_executable ( spex_demo_backslash         "Demo/spex_demo_backslash.c"         ${SPEX_DEMO_SOURCES} )
    add_executable ( spex_demo_cholesky_extended "Demo/spex_demo_cholesky_extended.c" ${SPEX_DEMO_SOURCES} )
    add_executable ( spex_demo_cholesky_simple   "Demo/spex_demo_cholesky_simple.c"   ${SPEX_DEMO_SOURCES} )
    add_executable ( spex_demo_lu_doub           "Demo/spex_demo_lu_doub.c"           ${SPEX_DEMO_SOURCES} )
    add_executable ( spex_demo_lu_extended       "Demo/spex_demo_lu_extended.c"       ${SPEX_DEMO_SOURCES} )
    add_executable ( spex_demo_lu_simple1        "Demo/spex_demo_lu_simple1.c"        ${SPEX_DEMO_SOURCES} )
    add_executable ( spex_demo_lu_simple2        "Demo/spex_demo_lu_simple2.c"        ${SPEX_DEMO_SOURCES} )
    add_executable ( spex_demo_Update_chol       "Demo/spex_demo_Update_chol.c"       ${SPEX_DEMO_SOURCES} )
    add_executable ( spex_demo_Update_lu         "Demo/spex_demo_Update_lu.c"         ${SPEX_DEMO_SOURCES} )

    # Libraries required for Demo programs
    target_link_libraries ( spex_demo_backslash         PUBLIC spex )
    target_link_libraries ( spex_demo_cholesky_extended PUBLIC spex )
    target_link_libraries ( spex_demo_cholesky_simple   PUBLIC spex )
    target_link_libraries ( spex_demo_lu_doub           PUBLIC spex )
    target_link_libraries ( spex_demo_lu_extended       PUBLIC spex )
    target_link_libraries ( spex_demo_lu_simple1        PUBLIC spex )
    target_link_libraries ( spex_demo_lu_simple2        PUBLIC spex )
    target_link_libraries ( spex_demo_Update_chol       PUBLIC spex )
    target_link_libraries ( spex_demo_Update_lu         PUBLIC spex )

else ( )

    message ( STATUS "Skipping the demos in SPEX/Demo" )

endif ( )

#-------------------------------------------------------------------------------
# report status
#-------------------------------------------------------------------------------

include ( SuiteSparseReport )

