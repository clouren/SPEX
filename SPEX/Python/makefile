#-------------------------------------------------------------------------------
# SPEX_Chol/Python/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX_Cholesky: (c) 2022, Chris Lourenco, United States Naval Academy,
# Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
# Texas A&M University. All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#-------------------------------------------------------------------------------

default: all

SUITESPARSE ?= $(realpath $(CURDIR)/../..)
include ../../SuiteSparse_config/SuiteSparse_config.mk

PATH_TO_SPEX=-Wl,-rpath,"../lib/"
CF =  $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(PATH_TO_SPEX) -O

I = -I../SPEX_Backslash/Include -I../../COLAMD/Include -I../../AMD/Include -I../../SuiteSparse_config -I../SPEX_Utilities/Source -I../SPEX_LU/Source -I../SPEX_Cholesky/Source -I../SPEX_Update/Source

LDLIBS += -lsuitesparseconfig -lcolamd -lamd \
           -lgmp -lm\
           -lmpfr -lm\
           -L../../lib -lspex -lm

all: SPEX_connect.o SPEX_connect.so

SPEX_connect.o: SPEX_connect.c
	$(CC) $(LDFLAGS) $(CF) $(I) -fPIC -c -o SPEX_connect.o SPEX_connect.c $(LDLIBS)
	
SPEX_connect.so: SPEX_connect.o
	$(CC) $(LDFLAGS) $(CF) $(I) -shared SPEX_connect.o -o SPEX_connect.so $(LDLIBS)
	


clean:
	- $(RM) SPEX_connect.so SPEX_connect.o

purge: distclean

distclean: clean
	- $(RM) -r SPEX_connect.so SPEX_connect.o *.a *.dSYM *.obj *.dll

