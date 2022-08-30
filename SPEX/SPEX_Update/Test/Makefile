#-------------------------------------------------------------------------------
# SPEX_Update/Test/Makefile: compile and run SPEX_update performance test 
#-------------------------------------------------------------------------------

# SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
# Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#-------------------------------------------------------------------------------

default: all

SUITESPARSE ?= $(realpath $(CURDIR)/../../..)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
CFLAGS += -Wall -Wextra -Wpedantic -Werror

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
I = -I../../../include -I../../SPEX_Utilities/Demo -I../../SPEX_Utilities/Source

LDLIBS += -lm -lgmp -lmpfr -lsuitesparseconfig -lspex -lglpk

all: lib ptest chol
	- ./ptest lp_afiro
	- ./cholupdate lp_afiro

runtest: lib ptest
	- ./ptest

runchol: lib chol
	- ./cholupdate

lib:
	( cd ../../build ; $(MAKE) )

ptest: lib ptest.c test.c test.h Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o ptest ptest.c test.c $(LDLIBS)

test_optimal: lib test_optimal.c test.c test.h Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o test_optimal test_optimal.c test.c $(LDLIBS)

openfile: lib openfile.c test.c test.h Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o openfile openfile.c test.c $(LDLIBS)

chol: lib cholupdate.c test.c test.h Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o cholupdate cholupdate.c test.c $(LDLIBS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r ptest *.a *.dSYM *.obj *.dll
