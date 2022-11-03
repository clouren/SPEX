#-------------------------------------------------------------------------------
# SPEX_LU/Demo/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX_LU: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
# Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#-------------------------------------------------------------------------------

default: all

SUITESPARSE ?= $(realpath $(CURDIR)/../../..)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
CFLAGS += -Wall -Wextra -Wpedantic -Werror

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
I = -I../../../include

LDLIBS += -lm -lsuitesparseconfig -lgmp -lmpfr -lcolamd -lamd -lspex

all: lib example example2 spexlu_demo spexlu_doub_demo
	- ./example
	- ./example2
	- ./spexlu_demo

lib:
	( cd ../../build ; $(MAKE) )

example: lib example.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example example.c $(LDLIBS)

example2: lib example2.c lu_demos.h lu_demos.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example2 example2.c lu_demos.c $(LDLIBS)

spexlu_demo: lib spexlu_demo.c lu_demos.h lu_demos.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o spexlu_demo spexlu_demo.c lu_demos.c $(LDLIBS)

spexlu_doub_demo: lib spexlu_doub_demo.c lu_demos.h lu_demos.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o spexlu_doub_demo spexlu_doub_demo.c lu_demos.c $(LDLIBS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r example example2 spexlu_demo spexlu_doub_demo *.a *.dSYM *.obj *.dll
