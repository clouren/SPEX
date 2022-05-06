#------------------------------------------------------------------------------
# SPEX_Left_LU/Makefile: run SPEX_Left_LU demo and tests
#------------------------------------------------------------------------------

# SPEX_Left_LU: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
# Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/../..)
export SUITESPARSE

# default: run the demos
default: demos

include ../../SuiteSparse_config/SuiteSparse_config.mk

# run the demos
demos:
	( cd Demo ; $(MAKE) )


# remove files not in the distribution, but keep compiled libraries
clean:
	( cd Demo ; $(MAKE) clean )
	- ( cd MATLAB/Source ; $(RM) *.o )
	- ( cd MATLAB    ; $(RM) *.o )
#	( cd Doc ; $(MAKE) clean )

# remove all files not in the distribution
purge:
	( cd Demo ; $(MAKE) purge )
	- ( cd MATLAB/Source ; $(RM) *.o *.mex* )
	- ( cd MATLAB    ; $(RM) *.o *.mex* )
#	( cd Doc ; $(MAKE) clean )

distclean: purge
