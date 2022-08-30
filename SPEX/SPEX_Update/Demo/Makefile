#-------------------------------------------------------------------------------
# SPEX_Update/Demo/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX_Update: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
# Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#-------------------------------------------------------------------------------

default: all

SUITESPARSE ?= $(realpath $(CURDIR)/../../..)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
CFLAGS += -Wall -Wextra -Wpedantic -Werror

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
I = -I../../../include -I../../SPEX_LU/Demo

LDLIBS += -lm -lsuitesparseconfig -lgmp -lmpfr -lcolamd -lamd -lspex

all: lib luupdate cholupdate
	- ./example_for_lu_update
	- ./example_for_chol_update

lib:
	( cd ../../build ; $(MAKE) )

luupdate: lib example_for_lu_update.c ../../SPEX_LU/Demo/demos.c ../../SPEX_LU/Demo/demos.h Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example_for_lu_update example_for_lu_update.c ../../SPEX_LU/Demo/demos.c $(LDLIBS)

cholupdate: lib example_for_chol_update.c ../../SPEX_LU/Demo/demos.c ../../SPEX_LU/Demo/demos.h Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example_for_chol_update example_for_chol_update.c ../../SPEX_LU/Demo/demos.c $(LDLIBS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r example_for_lu_update example_for_chol_update *.a *.dSYM *.obj *.dll
