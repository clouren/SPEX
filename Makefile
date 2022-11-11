#-------------------------------------------------------------------------------
# Makefile for all SPEX and its dependent packages
#-------------------------------------------------------------------------------

# edit this variable to pass options to cmake:
export CMAKE_OPTIONS ?=

# edit this variable to control parallel make:
export JOBS ?= 8

# do not modify this variable
export SUITESPARSE = $(CURDIR)

#-------------------------------------------------------------------------------

# Compile the default rules for each package.

# default: "make install" will install all libraries in /usr/local/lib
# and include files in /usr/local/include.  Not installed in SuiteSparse/lib.
default: library

# compile; "sudo make install" will install only in /usr/local
# (or whatever your CMAKE_INSTALL_PREFIX is)
library:
	( cd SuiteSparse_config && $(MAKE) )
	( cd AMD && $(MAKE) )
	( cd COLAMD && $(MAKE) )
	( cd SPEX && $(MAKE) )

# compile; "make install" only in  SuiteSparse/lib and SuiteSparse/include
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

# compile; "sudo make install" will install only in /usr/local
# (or whatever your CMAKE_INSTALL_PREFIX is)
both:
	( cd SuiteSparse_config && $(MAKE) both )
	( cd AMD && $(MAKE) both )
	( cd COLAMD && $(MAKE) both )
	( cd SPEX && $(MAKE) both )

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
	- $(RM) MATLAB_Tools/*/*.mex* MATLAB_Tools/*/*/*.mex*
	- $(RM) MATLAB_Tools/*/*.o    MATLAB_Tools/*/*/*.o
	- $(RM) -r include/* bin/* lib/*
	- ( cd SPEX && $(MAKE) purge )

# Remove all files not in the original distribution, but keep the libraries
clean:
	- ( cd SuiteSparse_config && $(MAKE) clean )
	- ( cd AMD && $(MAKE) clean )
	- ( cd COLAMD && $(MAKE) clean )
	- ( cd SPEX && $(MAKE) clean )

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

