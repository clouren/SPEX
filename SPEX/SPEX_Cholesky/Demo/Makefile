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
CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(PATH_TO_SPEX_UTIL) -O
I = -I../Include -I../Source -I../../../SuiteSparse_config -I../../../COLAMD/Include -I../../../AMD/Include -I../../SPEX_Util/Include -I../../SPEX_Util/Source -I../../SPEX_Util/Lib

LDLIBS += -lm -lsuitesparseconfig -lgmp -lmpfr -lcolamd -lamd -L../../SPEX_Util/Lib/ -lspexutil 
CS = ../Lib/libspexchol.a $(LDLIBS)


all: lib example_extended example_extended_debug example_simple #forCluster_extended
	- ./example_simple
	- ./example_extended
	#- ./forCluster_extended f ../ExampleMats/2.mat.txt s
	- echo "WARNING: Remember to remove the -g for test and releases in the Lib/Makefile"

lib:
	( cd ../Lib ; $(MAKE) )

example_extended: lib example_extended.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example_extended example_extended.c demos.c $(CS)

example_extended_debug: lib example_extended.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o example_extended_debug example_extended.c demos.c $(CS)

example_simple: lib example_simple.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o example_simple example_simple.c demos.c $(CS)

#forCluster_extended: lib forCluster_extended.c Makefile
#	$(CC) $(LDFLAGS) $(CF) $(I) -o forCluster_extended forCluster_extended.c demos.c $(CS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r example_simple example_extended example_extended_debug *.a *.dSYM *.obj *.dll
