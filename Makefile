#-------------------------------------------------------------------------------
# Makefile for SPEX and its dependent packages (AMD, COLAMD, SuiteSparse_config)
#-------------------------------------------------------------------------------

# Copyright (c) 2023-2024, Timothy A. Davis, All Rights Reserved. FIXME
# Just this particular file is under the Apache-2.0 license; each package has
# its own license.
# SPDX-License-Identifier: Apache-2.0

# edit this variable to pass options to cmake:
export CMAKE_OPTIONS ?=

# edit this variable to control parallel make:
export JOBS ?= 8

# do not modify this variable
export SUITESPARSE = $(CURDIR)

#-------------------------------------------------------------------------------

# Compile the default rules for each package.

# default: compile and install in SPEX/lib and SPEX/include
default: local install

# compile; "sudo make install" will install only in /usr/local
# (or whatever your CMAKE_INSTALL_PREFIX is)
library:
	( cd SuiteSparse_config && $(MAKE) )
	( cd AMD && $(MAKE) )
	( cd COLAMD && $(MAKE) )
	( cd SPEX && $(MAKE) )

# compile; "make install" only in SPEX/lib and SPEX/include
local:
	( cd SuiteSparse_config && $(MAKE) local )
	( cd AMD && $(MAKE) local )
	( cd COLAMD && $(MAKE) local )
	( cd SPEX && $(MAKE) local )

# compile; "sudo make install" will install only in /usr/local
# (or whatever your CMAKE_INSTALL_PREFIX is)
global:
	( cd SuiteSparse_config && $(MAKE) global )
	( cd AMD && $(MAKE) global )
	( cd COLAMD && $(MAKE) global )
	( cd SPEX && $(MAKE) global )

# install all packages.  Location depends on prior "make", "make global" etc
install:
	( cd SuiteSparse_config && $(MAKE) install )
	( cd AMD && $(MAKE) install )
	( cd COLAMD && $(MAKE) install )
	( cd SPEX && $(MAKE) install )

# uninstall all packages
uninstall:
	( cd SuiteSparse_config && $(MAKE) uninstall )
	( cd AMD && $(MAKE) uninstall )
	( cd COLAMD && $(MAKE) uninstall )
	( cd SPEX && $(MAKE) uninstall )

# Remove all files not in the original distribution
distclean: purge

# Remove all files not in the original distribution
purge:
	- ( cd SuiteSparse_config && $(MAKE) purge )
	- ( cd AMD && $(MAKE) purge )
	- ( cd COLAMD && $(MAKE) purge )
	- ( cd SPEX && $(MAKE) purge )
	- $(RM) -r include/* bin/* lib/*

clean: purge

# Run all demos
demos:
	- ( cd SuiteSparse_config && $(MAKE) demos )
	- ( cd AMD && $(MAKE) demos )
	- ( cd COLAMD && $(MAKE) demos )
	- ( cd SPEX && $(MAKE) demos )

# Create the PDF documentation
docs:
	( cd AMD && $(MAKE) docs )
	( cd SPEX && $(MAKE) docs )

# statement coverage (Linux only); this requires a lot of time.
cov: local install
	( cd SPEX && $(MAKE) cov )

debug:
	( cd SuiteSparse_config && $(MAKE) debug )
	( cd AMD && $(MAKE) debug )
	( cd COLAMD && $(MAKE) debug )
	( cd SPEX && $(MAKE) debug )

