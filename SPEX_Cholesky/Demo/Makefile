#-------------------------------------------------------------------------------
# SPEX_Chol/Demo/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
# Erick Moreno-Centeno, Timothy A. Davis, Texas A&M University.  All Rights Reserved.  
# See SPEX_Cholesky/License for the license.

#-------------------------------------------------------------------------------

default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
#CFLAGS += -Wall -Wextra -Wpedantic -Werror

PATH_TO_SPEX_UTIL=-Wl,-rpath,"../../SPEX_Util/Lib/"
CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(PATH_TO_SPEX_UTIL) -O
I = -I../Include -I../Source -I../../SuiteSparse_config -I../../COLAMD/Include -I../../AMD/Include -I../../SPEX_Util/Include -I../../SPEX_Util/Source -I../../SPEX_Util/Lib

LDLIBS += -lm -lsuitesparseconfig -lgmp -lmpfr -lcolamd -lamd -L../../SPEX_Util/Lib/ -lspexutil 
CS = ../Lib/libspexchol.a $(LDLIBS)


all: lib SPEX_Chol SPEX_Chol_debug
	- ./SPEX_Chol

lib:
	( cd ../Lib ; $(MAKE) )

SPEX_Chol: lib SPEX_Chol.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o SPEX_Chol SPEX_Chol.c $(CS)

SPEX_Chol_debug: lib SPEX_Chol.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o SPEX_Chol_debug SPEX_Chol.c $(CS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r SPEX_Chol SPEX_Chol_debug *.a *.dSYM *.obj *.dll
