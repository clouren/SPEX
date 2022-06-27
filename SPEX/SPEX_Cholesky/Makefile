#------------------------------------------------------------------------------
# SPEX_Chol/Makefile: run SPEX Cholesky demo and tests
#------------------------------------------------------------------------------

# SPEX_Cholesky: (c) 2022, Chris Lourenco, United States Naval Academy, 
# Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
# Texas A&M University. All Rights Reserved. 
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
	#( cd Tcov ; $(MAKE) clean )
	- ( cd MATLAB/Source ; $(RM) *.o )
	- ( cd MATLAB    ; $(RM) *.o )
#	( cd Doc ; $(MAKE) clean )

# remove all files not in the distribution
purge:
	( cd Demo ; $(MAKE) purge )
	#( cd Tcov ; $(MAKE) purge )
	- ( cd MATLAB/Source ; $(RM) *.o *.mex* )
	- ( cd MATLAB    ; $(RM) *.o *.mex* )
#	( cd Doc ; $(MAKE) clean )

distclean: purge
