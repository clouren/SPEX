#-------------------------------------------------------------------------------
# SPEX_Chol/Demo/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX_Cholesky: (c) 2022, Chris Lourenco, United States Naval Academy,
# Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
# Texas A&M University. All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#-------------------------------------------------------------------------------

default: all

SUITESPARSE ?= $(realpath $(CURDIR)/../../..)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
CFLAGS += -Wall -Wextra -Wpedantic -Werror

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
I = -I../../../include -I../../SPEX_Utilities/Source

LDLIBS += -lm -lsuitesparseconfig -lgmp -lmpfr -lcolamd -lamd -lspex


all: lib example_extended example_extended_debug example_simple
	- ./example_simple
	- ./example_extended

lib:
	( cd ../../build ; $(MAKE) )

example_extended: lib example_extended.c chol_demos.h chol_demos.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example_extended example_extended.c chol_demos.c $(LDLIBS)

example_extended_debug: lib example_extended.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example_extended_debug example_extended.c chol_demos.c $(LDLIBS)

example_simple: lib example_simple.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example_simple example_simple.c chol_demos.c $(LDLIBS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r example_simple example_extended example_extended_debug *.a *.dSYM *.obj *.dll
