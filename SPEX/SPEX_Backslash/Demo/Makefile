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
CFLAGS += -Wall -Wextra -Wpedantic -Werror

CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O
I = -I../../../include

LDLIBS += -lm -lsuitesparseconfig -lgmp -lmpfr -lcolamd -lamd -lspex


# FIXME: the 2nd command reports an error when openning the RHS file,
# then generates its own RHS.
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
	- $(RM) -r *.a *.dSYM *.obj *.dll
