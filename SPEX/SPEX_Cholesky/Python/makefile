#-------------------------------------------------------------------------------
# SPEX_Chol/Demo/Makefile: compile and run the demos
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

PATH_TO_SPEX_UTIL=-Wl,-rpath,"../../SPEX_Util/Lib/"
PATH_TO_SPEX_CHOL=-Wl,-rpath,"../SPEX/SPEX_Cholesky/Lib/"
CF =  $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(PATH_TO_SPEX_UTIL) $(PATH_TO_SPEX_CHOL) -O
I = -I../Include -I../Source -I../Lib -I../../../SuiteSparse_config -I../../../COLAMD/Include -I../../../AMD/Include -I../../SPEX_Util/Include -I../../SPEX_Util/Source -I../../SPEX_Util/Lib

LDLIBS += -lsuitesparseconfig -lcolamd -lamd \
          -lgmp -lmpfr \
           -L../../SPEX_Util/Lib/ -lspexutil -lm \
           -L../Lib -lspexchol \
           -lm -lsuitesparseconfig -lcolamd -lamd -lgmp -lmpfr 
#CS = ../Lib/libspexchol.a $(LDLIBS)




all: demo demo_debug example.o libraryE.so
	- ./demo
	- ./demo_debug


demo: example.c 
	$(CC) $(LDFLAGS) $(CF) $(I) -o demo main.c example.c $(LDLIBS)
	
demo_debug: example.c 
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o demo_debug main.c example.c $(LDLIBS)

example.o: example.c
	$(CC) $(LDFLAGS) $(CF) $(I) -fPIC -c -o example.o example.c $(LDLIBS)
	
libraryE.so: example.o 
	$(CC) $(LDFLAGS) $(CF) $(I) -shared example.o -o libraryE.so $(LDLIBS)
	


clean:
	- $(RM) demo example libraryE.so example.o

purge: distclean

distclean: clean
	- $(RM) -r SPEX_Chol SPEX_Chol_debug *.a *.dSYM *.obj *.dll

