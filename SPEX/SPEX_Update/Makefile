#------------------------------------------------------------------------------
# SPEX_LU_Update/Makefile: compile, install, and run SPEX_LU_Update demo
# and tests
#------------------------------------------------------------------------------

# SPEX_LU_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
# Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
# SPEX_LU_Update/License for the license.

#------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/../..)
#export SUITESPARSE

# default: compile the dynamic library, and run the demos
default: C

include ../../SuiteSparse_config/SuiteSparse_config.mk

# compile the dynamic library, and then run the demos and statement coverage
demos: all

# compile the dynamic library only, and then run the demos
C:
	( cd Lib ; $(MAKE) )
	#( cd Demo ; $(MAKE) )

# compile the dynamic library, and then run the demos and statement coverage
all: C cov

# compile the dynamic library only
library:
	( cd Lib ; $(MAKE) )

# compile the static library only
static:
	( cd Lib ; $(MAKE) static )

# statement coverage test
cov:
	( cd Tcov ; $(MAKE) )

# remove files not in the distribution, but keep compiled libraries
clean:
	( cd Lib ; $(MAKE) clean )
	( cd Demo ; $(MAKE) clean )
	( cd Tcov ; $(MAKE) clean )
	#- ( cd MATLAB/Source ; $(RM) *.o )
	#- ( cd MATLAB    ; $(RM) *.o )
	#( cd Doc ; $(MAKE) clean )

# remove all files not in the distribution
purge:
	( cd Lib ; $(MAKE) purge )
	( cd Tcov ; $(MAKE) purge )
	( cd Demo ; $(MAKE) purge )
	#- ( cd MATLAB/Source ; $(RM) *.o *.mex* )
	#- ( cd MATLAB    ; $(RM) *.o *.mex* )
	#( cd Doc ; $(MAKE) clean )

distclean: purge

# install the library
install: library
	( cd Lib ; $(MAKE) install )

# uninstall the library
uninstall:
	( cd Lib ; $(MAKE) uninstall )

