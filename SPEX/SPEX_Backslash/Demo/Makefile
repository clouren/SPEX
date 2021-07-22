#-------------------------------------------------------------------------------
# SPEX/SPEX_Backslash/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX: (c) 2021
# This is a general purpose makefile that contains all of the components of SPEX
#-------------------------------------------------------------------------------

default: all

# Relative path to Suitesparse
SUITESPARSE ?= $(realpath $(CURDIR)/../../../)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
#CFLAGS += -Wall -Wextra -Wpedantic -Werror

# Here is the path to all SPEX routines. Any new SPEX Routine can be added here
PATH_TO_SPEX_UTIL=-Wl,-rpath,"../../SPEX_Util/Lib/"
PATH_TO_SPEX_CHOL=-Wl,-rpath,"../../SPEX_Cholesky/Lib/"
PATH_TO_SPEX_LEFT=-Wl,-rpath,"../../SPEX_Left_LU/Lib/"
PATH_TO_SPEX_BACK=-Wl,-rpath,"../Lib/"

# CF Needs to include the CFLAGS plus the path to all SPEX routines
CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(PATH_TO_SPEX_UTIL) $(PATH_TO_SPEX_CHOL) \
     $(PATH_TO_SPEX_LEFT) $(PATH_TO_SPEX_BACK) -O

# Here we include all the dependencies. The first line is SuiteSparse COLAD and AMD
# The second line is SPEX Util 
# The third line is SPEX Left LU
# The fourth line is SPEX Cholesky 
# The final line is SPEX_Backslash
I = -I../../../SuiteSparse_config -I../../../COLAMD/Include -I../../../AMD/Include \
    -I../../SPEX_Util/Include -I../../SPEX_Util/Source -I../../SPEX_Util/Lib \
    -I../../SPEX_Left_LU/Include -I../../SPEX_Left_LU/Source -I../../SPEX_Left_LU/Lib \
    -I../../SPEX_Cholesky/Include -I../../SPEX_Cholesky/Source -I../../SPEX_Cholesky/Lib \
    -I../Include -I../Source -I../Lib

# Include all the lib files. This includes the suitesparse, gmp, and spex libs
# First line is Suitesparse 
# The second line is gmp and mpfr
# The third line is SPEX Util 
# The fourth line is SPEX Left LU
# The fifth line is SPEX Cholesky 
# The sixth line is SPEX WAMF 

#INSERT SPEX_BACKSLASH WHEN MADE!

# The last line is m
LDLIBS += -lsuitesparseconfig -lcolamd -lamd \
          -lgmp -lmpfr \
          -L../../SPEX_Util/Lib/ -lspexutil -lm \
          -L../../SPEX_Left_LU/Lib/ -lspexleftlu -lm \
          -L../../SPEX_Cholesky/Lib/ -lspexchol -lm \
          -L../Lib/ -lspexbackslash -lm \
          -lm -lsuitesparseconfig -lcolamd -lamd -lgmp -lmpfr 

#CS = ../Lib/libspexwamf.a $(LDLIBS)

all: SPEX_Backslash_Demo
	- ./SPEX_Backslash_Demo
	- ./SPEX_Backslash_Demo f ../ExampleMats/2209.mat.txt s

lib:
	( cd ../ ; $(MAKE lib) )

SPEX_Backslash_Demo: lib SPEX_Backslash_Demo.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o SPEX_Backslash_Demo SPEX_Backslash_Demo.c demos.c $(LDLIBS)


clean:
	- $(RM) SPEX_Backslash_Demo

purge: distclean

distclean: clean
	- $(RM) -r SPEX_Chol SPEX_Chol_debug *.a *.dSYM *.obj *.dll
