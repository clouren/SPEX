#-------------------------------------------------------------------------------
# SPEX_QR/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
# Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
# All Rights Reserved.  See SPEX_Cholesky/License for the license.

#-------------------------------------------------------------------------------

default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
#CFLAGS += -Wall -Wextra -Wpedantic -Werror

PATH_TO_SPEX_UTIL=-Wl,-rpath,"../SPEX_Util/Lib/"
CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(PATH_TO_SPEX_UTIL) -O2
CF_DEBUG = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(PATH_TO_SPEX_UTIL) -g
I = -I../../SuiteSparse_config -I../../COLAMD/Include -I../../AMD/Include -I../SPEX_Util/Include -I../SPEX_Util/Source -I../SPEX_Util/Lib -I./

LDLIBS += -lm  -lcolamd -lamd -L../SPEX_Util/Lib/ -lspexutil -lgmp -lmpfr -lsuitesparseconfig
CS = $(LDLIBS)


all: SPEX_QR_dense SPEX_QR_dense_debug

SPEX_QR_dense: SPEX_QR_dense.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o SPEX_QR_dense SPEX_QR_dense.c $(CS)

SPEX_QR_dense_debug: SPEX_QR_dense.c Makefile
	$(CC) $(LDFLAGS) $(CF_DEBUG) $(I) -o debug_SPEX_QR_dense SPEX_QR_dense.c $(CS)

clean:
	- $(RM) *.o SPEX_QR_dense

purge: distclean

distclean: clean
	- $(RM) -r SPEX_Chol SPEX_Chol_debug *.a *.dSYM *.obj *.dll
