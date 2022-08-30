#------------------------------------------------------------------------------
# SPEX_Backslash/Makefile: compile, install, and run SPEX_Backslash demo and tests
#------------------------------------------------------------------------------

# SPEX_Backslash: (c) 2021, Chris Lourenco, United States Naval Academy, 
# Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
# Texas A&M University. All Rights Reserved. 
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#------------------------------------------------------------------------------

SUITESPARSE ?= $(realpath $(CURDIR)/../..)
export SUITESPARSE

# default: compile the dynamic library, and run the demos
default: C

include ../../SuiteSparse_config/SuiteSparse_config.mk

# compile the dynamic library, and then run the demos and statement coverage
demos: all

# compile the dynamic library only, and then run the demos
C:
	( cd Demo ; $(MAKE) )

# compile the dynamic library, and then run the demos and statement coverage
all: C cov

# remove files not in the distribution, but keep compiled libraries
clean:
	( cd Demo ; $(MAKE) clean )

# remove all files not in the distribution
purge:
	( cd Demo ; $(MAKE) purge )

distclean: purge


