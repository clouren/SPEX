#-------------------------------------------------------------------------------
# SuiteSparse/SPEX/Makefile: Makefile for SPEX
#-------------------------------------------------------------------------------

# SPEX: (c) 2019-2022, Chris Lourenco, Jinhao Chen,
# Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
# All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#-------------------------------------------------------------------------------

# A simple Makefile for SPEX, which relies on cmake to do the
# actual build.  All the work is done in cmake so this Makefile is just for
# convenience.

# To compile with an alternate compiler:
#
#       make CC=gcc CXX=g++
#
# To compile/install for system-wide usage:
#
#       make
#       sudo make install
#
# To compile/install for local usage (SuiteSparse/lib and SuiteSparse/include):
#
#       make local
#       make install
#
# To clean up the files:
#
#       make clean

JOBS ?= 8

default: library

# default is to install only in /usr/local
library:
	( cd build && cmake $(CMAKE_OPTIONS) .. && cmake --build . -j${JOBS} )

# install only in SuiteSparse/lib and SuiteSparse/include
local:
	( cd build && cmake $(CMAKE_OPTIONS) -DGLOBAL_INSTALL=0 -DLOCAL_INSTALL=1 .. && cmake --build . -j${JOBS} )

# install only in /usr/local (default)
global:
	( cd build && cmake $(CMAKE_OPTIONS) -DGLOBAL_INSTALL=1 -DLOCAL_INSTALL=0 .. && cmake --build . -j${JOBS} )

# install in SuiteSparse/lib both and SuiteSparse/include and /usr/local
both:
	( cd build && cmake $(CMAKE_OPTIONS) -DGLOBAL_INSTALL=1 -DLOCAL_INSTALL=1 .. && cmake --build . -j${JOBS} )

debug:
	( cd build && cmake $(CMAKE_OPTIONS) -DCMAKE_BUILD_TYPE=Debug .. && cmake --build . -j${JOBS} )

all: library

# FIXME: "s" is a bad idea to generate an all 1's RHS.
# FIXME: what is 2.mat.soln.txt?  It's being used as a RHS input, not a solution
demos:
	( cd build && cmake $(CMAKE_OPTIONS) -DDEMO=1 .. && cmake --build . -j${JOBS} )
	./build/spex_demo_lu_simple1
	./build/spex_demo_lu_simple2          ExampleMats/10teams_mat.txt  ExampleMats/10teams_v.txt
	./build/spex_demo_lu_extended       f ExampleMats/10teams_mat.txt  ExampleMats/10teams_v.txt
	./build/spex_demo_lu_doub           f ExampleMats/10teams_mat.txt  ExampleMats/10teams_v.txt
	./build/spex_demo_backslash_demo    f ExampleMats/10teams_mat.txt  ExampleMats/10teams_v.txt
	./build/spex_demo_backslash_demo    f ExampleMats/2209.mat.txt     s
	./build/spex_demo_cholesky_simple   f ExampleMats/2.mat.txt        ExampleMats/2.mat.soln.txt
	./build/spex_demo_cholesky_extended f ExampleMats/2.mat.txt        ExampleMats/2.mat.soln.txt
	./build/spex_demo_Update_lu         f ExampleMats/10teams_mat.txt  ExampleMats/10teams_v.txt
	./build/spex_demo_Update_chol       f ExampleMats/872_mat.txt      s

# FIXME remove olddemos target
olddemos:
	( cd SPEX_LU/Demo && ../../build/example  )
	( cd SPEX_LU/Demo && ../../build/example2 )
	( cd SPEX_LU/Demo && ../../build/spexlu_demo )
	( cd SPEX_LU/Demo && ../../build/spexlu_doub_demo )
	( cd SPEX_Backslash/Demo && ../../build/SPEX_Backslash_Demo )
	( cd SPEX_Backslash/Demo && ../../build/SPEX_Backslash_Demo f ../../ExampleMats/2209.mat.txt s )
	( cd SPEX_Cholesky/Demo && ../../build/example_simple )
	( cd SPEX_Cholesky/Demo && ../../build/example_extended )
	( cd SPEX_Update/Demo && ../../build/example_for_lu_update )
	( cd SPEX_Update/Demo && ../../build/example_for_chol_update )

cov:
	( cd Tcov && $(MAKE) )

# just compile after running cmake; do not run cmake again
remake:
	( cd build && cmake --build . -j${JOBS} )

# just run cmake to set things up
setup:
	( cd build && cmake $(CMAKE_OPTIONS) .. )

install:
	( cd build && cmake --install . )

# remove any installed libraries and #include files
uninstall:
	- xargs rm < build/install_manifest.txt

# remove all files not in the distribution
clean:
	- $(RM) -rf build/* Config/*.tmp MATLAB/*.o MATLAB/*.mex*
	( cd Tcov && $(MAKE) clean )

purge: clean

distclean: clean

docs:
	( cd Doc && $(MAKE) )
