#-------------------------------------------------------------------------------
# SPEX_LU/Demo/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX_LU: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
# Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
# SPEX_LU/License for the license.

#-------------------------------------------------------------------------------

default: all

include ../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
CFLAGS += -Wall -Wextra -Wpedantic -Werror

# CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
I = -I../Include -I../Source -I../../SuiteSparse_config -I../../COLAMD/Include -I../../AMD/Include

# LDFLAGS = -L../../lib

LDLIBS += -lm -lgmp -lmpfr -lcolamd -lamd -lspexlu
CS = $(LDLIBS)


all: lib example 
	- ./example

lib:
	( cd ../Lib ; $(MAKE) )

example: lib example.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -o example example.c $(CS)
clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r example example2 example3 example4 example5 SPEXLU *.a *.dSYM *.obj *.dll
