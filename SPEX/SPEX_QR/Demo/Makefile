#-------------------------------------------------------------------------------
# SPEX_QR/Demo/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX_QR: (c) 2021, Chris Lourenco, US Naval Academy, All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#-------------------------------------------------------------------------------

default: all

SUITESPARSE ?= $(realpath $(CURDIR)/../../..)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
#CFLAGS += -Wall -Wextra -Wpedantic -Werror

PATH_TO_SPEX_UTIL=-Wl,-rpath,"../../SPEX_Util/Lib/"
CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(PATH_TO_SPEX_UTIL) -O
I = -I../Include -I../Source -I../../../SuiteSparse_config -I../../../COLAMD/Include -I../../../AMD/Include -I../../SPEX_Util/Include -I../../SPEX_Util/Source -I../../SPEX_Util/Lib

LDLIBS += -lm -lsuitesparseconfig -lgmp -lmpfr -lcolamd -lamd -L../../SPEX_Util/Lib/ -lspexutil 
CS = ../Lib/libspexqr.a $(LDLIBS)


all: lib SPEX_QR_dense_demo SPEX_QR_dense_demo_debug
	- ./SPEX_QR_dense 100 100 100 10 100

lib:
	( cd ../Lib ; $(MAKE) )

SPEX_QR_dense_demo: lib SPEX_QR_dense_demo.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o SPEX_QR_dense SPEX_QR_dense_demo.c $(CS)

SPEX_QR_dense_demo_debug: lib SPEX_QR_dense_demo.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o SPEX_QR_dense_debug SPEX_QR_dense_demo.c $(CS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r SPEX_QR_dense SPEX_QR_dense_debug *.a *.dSYM *.obj *.dll
