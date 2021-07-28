#-------------------------------------------------------------------------------
# SPEX_LU_Update/Test/Makefile: compile and run SPEX_LU_update performance test 
#-------------------------------------------------------------------------------

# SPEX_LU_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
# Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
# SPEX_LU_Update/License for the license.

#-------------------------------------------------------------------------------

default: all

SUITESPARSE ?= $(realpath $(CURDIR)/../../..)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
CFLAGS += -Wall -Wextra -Wpedantic -Werror

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
I = -I../Include -I../Lib -I../../../SuiteSparse_config -I../../SPEX_Util/Include -I../../SPEX_Util/Lib -I../../SPEX_Left_LU/Include -I../../SPEX_Left_LU/Lib -I../../SPEX_Cholesky/Include -I../../SPEX_Cholesky/Lib -I../../../lib -I/opt/ibm/ILOG/CPLEX_Studio201/cplex/include -L/opt/ibm/ILOG/CPLEX_Studio201/cplex/bin/x86-64_linux/
#-L/opt/ibm/ILOG/CPLEX_Studio201/cplex/lib/x86-64_linux/static_pic

LDLIBS += -lm -lgmp -lsuitesparseconfig -lspexutil -lspexupdate -lspexleftlu -lspexchol -lglpk -lcplex2010

all: lib ptest
	- ./ptest

runchol: lib chol
	- ./cholupdate

lib:
	( cd ../Lib ; $(MAKE) )

cptest: lib cplex_ptest.c cplex_test.c cplex_test.h Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o cplex_ptest cplex_ptest.c cplex_test.c $(LDLIBS)

ptest: lib ptest.c test.c test.h Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o ptest ptest.c test.c $(LDLIBS)

chol: lib cholupdate.c test.c test.h Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o cholupdate cholupdate.c test.c $(LDLIBS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r ptest *.a *.dSYM *.obj *.dll
