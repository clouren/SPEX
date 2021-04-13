#-------------------------------------------------------------------------------
# SPEX/Tests/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX: (c) 2021
# This is a general purpose makefile that contains all of the components of SPEX
#-------------------------------------------------------------------------------

default: all

# Relative path to Suitesparse
SUITESPARSE ?= $(realpath $(CURDIR)/../)
include ../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
#CFLAGS += -Wall -Wextra -Wpedantic -Werror

# Here is the path to all SPEX routines. Any new SPEX Routine can be added here
PATH_TO_SPEX_UTIL=-Wl,-rpath,"../SPEX/SPEX_Util/Lib/"
PATH_TO_SPEX_CHOL=-Wl,-rpath,"../SPEX/SPEX_Cholesky/Lib/"
PATH_TO_SPEX_LEFT=-Wl,-rpath,"../SPEX/SPEX_Left_LU/Lib/"
PATH_TO_SPEX_WAMF=-Wl,-rpath,"../SPEX/SPEX_WAMF/Lib/"

# CF Needs to include the CFLAGS plus the path to all SPEX routines
CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(PATH_TO_SPEX_UTIL) $(PATH_TO_SPEX_CHOL) \
     $(PATH_TO_SPEX_LEFT) $(PATH_TO_SPEX_WAMF) -O

# Here we include all the dependencies. The first line is SuiteSparse COLAD and AMD
# The second line is SPEX Util 
# The third line is SPEX Left LU
# The fourth line is SPEX Cholesky 
# The final line is SPEX WAMF 
I = -I../SuiteSparse_config -I../COLAMD/Include -I../AMD/Include \
    -I../SPEX/SPEX_Util/Include -I../SPEX/SPEX_Util/Source -I../SPEX/SPEX_Util/Lib \
    -I../SPEX/SPEX_Left_LU/Include -I../SPEX/SPEX_Left_LU/Source -I../SPEX/SPEX_Left_LU/Lib \
    -I../SPEX/SPEX_Cholesky/Include -I../SPEX/SPEX_Cholesky/Source -I../SPEX/SPEX_Cholesky/Lib \
    -I../SPEX/SPEX_WAMF/Include -I../SPEX/SPEX_WAMF/Source -I../SPEX/SPEX_WAMF/Lib

# Include all the lib files. This includes the suitesparse, gmp, and spex libs
# First line is Suitesparse 
# The second line is gmp and mpfr
# The third line is SPEX Util 
# The fourth line is SPEX Left LU
# The fifth line is SPEX Cholesky 
# The sixth line is SPEX WAMF 
# The last line is m
LDLIBS += -lsuitesparseconfig -lcolamd -lamd \
          -lgmp -lmpfr \
          -L../SPEX/SPEX_Util/Lib/ -lspexutil -lm \
          -L../SPEX/SPEX_Left_LU/Lib/ -lspexleftlu -lm \
          -L../SPEX/SPEX_Cholesky/Lib/ -lspexchol -lm \
          -L../SPEX/SPEX_WAMF/Lib/ -lspexwamf \
          -lm -lsuitesparseconfig -lcolamd -lamd -lgmp -lmpfr 

#CS = ../Lib/libspexwamf.a $(LDLIBS)

all: SPEX_WAMF_Unsymmetric_Basis SPEX_WAMF_Unsymmetric_SJ SPEX_WAMF_SPD
	- ./SPEX_WAMF_Unsymmetric_Basis
	- ./SPEX_WAMF_Unsymmetric_SJ
	- ./SPEX_WAMF_SPD

lib:
	( cd ../ ; $(MAKE lib) )

SPEX_WAMF_Unsymmetric_Basis: lib SPEX_WAMF_Unsymmetric_Basis.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o SPEX_WAMF_Unsymmetric_Basis SPEX_WAMF_Unsymmetric_Basis.c demos.c $(LDLIBS)


SPEX_WAMF_Unsymmetric_SJ: lib SPEX_WAMF_Unsymmetric_SJ.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o SPEX_WAMF_Unsymmetric_SJ SPEX_WAMF_Unsymmetric_SJ.c demos.c $(LDLIBS)

SPEX_WAMF_SPD: lib SPEX_WAMF_SPD.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o SPEX_WAMF_SPD SPEX_WAMF_SPD.c demos.c $(LDLIBS)

clean:
	- $(RM) *SPEX_WAMF_Unsymmetric_Basis SPEX_WAMF_Unsymmetric_SJ SPEX_WAMF_SPD

purge: distclean

distclean: clean
	- $(RM) -r SPEX_Chol SPEX_Chol_debug *.a *.dSYM *.obj *.dll
