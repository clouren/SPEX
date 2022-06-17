#-------------------------------------------------------------------------------
# SPEX_Chol/Python/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
# Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
# All Rights Reserved.  See SPEX_Cholesky/License for the license.

#-------------------------------------------------------------------------------

default: all

SUITESPARSE ?= $(realpath $(CURDIR)/../../..)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
#CFLAGS += -Wall -Wextra -Wpedantic -Werror

#PATH_TO_SPEX_UTIL=-Wl,-rpath,"../../SPEX_Util/Lib/"
#PATH_TO_SPEX_CHOL=-Wl,-rpath,"../SPEX/SPEX_Cholesky/Lib/"
PATH_TO_SPEX=-Wl,-rpath,"../../lib/"
CF =  $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(PATH_TO_SPEX) -O #$(PATH_TO_SPEX_UTIL) $(PATH_TO_SPEX_CHOL) -O
I = -I../Include -I../Source -I../Lib -I../../../SuiteSparse_config -I../../../COLAMD/Include -I../../../AMD/Include -I../../SPEX_Util/Include -I../../SPEX_Util/Source -I../../SPEX_Util/Lib -I../../lib -I../../SPEX_Backslash/Include

LDLIBS += -lsuitesparseconfig -lcolamd -lamd \
          -lgmp -lmpfr \
           -L../../lib -lspex -lm \
           -L../Lib -lspexchol \
           -lm -lsuitesparseconfig -lcolamd -lamd -lgmp -lmpfr 
#CS = ../Lib/libspexchol.a $(LDLIBS)




all: demo demo_debug SPEX_Chol_connect.o SPEX_Chol_connect.so
	- ./demo
	- ./demo_debug


demo: SPEX_Chol_connect.c 
	$(CC) $(LDFLAGS) $(CF) $(I) -o demo main.c SPEX_Chol_connect.c $(LDLIBS)
	
demo_debug: SPEX_Chol_connect.c 
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o demo_debug main.c SPEX_Chol_connect.c $(LDLIBS)

SPEX_Chol_connect.o: SPEX_Chol_connect.c
	$(CC) $(LDFLAGS) $(CF) $(I) -fPIC -c -o SPEX_Chol_connect.o SPEX_Chol_connect.c $(LDLIBS)
	
SPEX_Chol_connect.so: SPEX_Chol_connect.o 
	$(CC) $(LDFLAGS) $(CF) $(I) -shared SPEX_Chol_connect.o -o SPEX_Chol_connect.so $(LDLIBS)
	


clean:
	- $(RM) demo demo_debug SPEX_Chol_connect.so SPEX_Chol_connect.o

purge: distclean

distclean: clean
	- $(RM) -r example example_debug *.a *.dSYM *.obj *.dll

