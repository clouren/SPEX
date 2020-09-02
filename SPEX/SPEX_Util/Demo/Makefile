#-------------------------------------------------------------------------------
# SPEX_LU/Demo/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
# Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
# SPEX_LU/License for the license.

#-------------------------------------------------------------------------------

default: all

include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
CFLAGS += -Wall -Wextra -Wpedantic -Werror

# CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
I = -I../Include -I../../../include -I../Source -I../../../SuiteSparse_config -I../../../COLAMD/Include -I../../../AMD/Include -I../Lib

# LDFLAGS = -L../../lib

LDLIBS += -lm -lsuitesparseconfig -lgmp -lmpfr -lcolamd -lamd -lspexutil
CS = $(LDLIBS)


all: lib demo
	- ./demo

lib:
	( cd ../Lib ; $(MAKE) )

demo: lib demo.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o demo demo.c $(CS)
clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r demo
