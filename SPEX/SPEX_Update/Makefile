#------------------------------------------------------------------------------
# SPEX_Update/Makefile: run SPEX_Update demo and tests
#------------------------------------------------------------------------------

# SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
# Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/../..)
#export SUITESPARSE

# default: run the demos
default: demos

include ../../SuiteSparse_config/SuiteSparse_config.mk

# run the demos
demos:
	( cd Demo ; $(MAKE) )

# remove files not in the distribution, but keep compiled libraries
clean:
	( cd Demo ; $(MAKE) clean )
	#( cd Test ; $(MAKE) clean )
	#- ( cd MATLAB/Source ; $(RM) *.o )
	#- ( cd MATLAB    ; $(RM) *.o )

# remove all files not in the distribution
purge:
	( cd Demo ; $(MAKE) purge )
	#( cd Test ; $(MAKE) purge )
	#- ( cd MATLAB/Source ; $(RM) *.o *.mex* )
	#- ( cd MATLAB    ; $(RM) *.o *.mex* )

distclean: purge
