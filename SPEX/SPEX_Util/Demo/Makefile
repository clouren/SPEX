#-------------------------------------------------------------------------------
# SPEX_Util/Demo/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX_Util: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
# Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#-------------------------------------------------------------------------------

default: all

SUITESPARSE ?= $(realpath $(CURDIR)/../../..)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
CFLAGS += -Wall -Wextra -Wpedantic -Werror

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
I = -I../../../include -I../Source/

LDLIBS += -lm -lsuitesparseconfig -lgmp -lmpfr -lcolamd -lamd -lspex

all: lib example
	- ./example

lib:
	( cd ../../build ; $(MAKE) )

example: lib example.c my_safe_gmp.h my_safe_gmp.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example example.c my_safe_gmp.c $(LDLIBS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r example *.a *.dSYM *.obj *.dll
