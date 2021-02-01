#-------------------------------------------------------------------------------
# SPEX_WAMF/Demo/Makefile: compile and run the demos
#-------------------------------------------------------------------------------

# SPEX_WAMF: (c) 2021, Chris Lourenco US Naval Academy, Erick Moreno-Centeno, Texas
# A&M University. All Rights Reserved.  See SPEX_WAMF/License for the license.

#-------------------------------------------------------------------------------

default: all

include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
#CFLAGS += -Wall -Wextra -Wpedantic -Werror

PATH_TO_SPEX_UTIL=-Wl,-rpath,"../../SPEX_Util/Lib/"
CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) $(PATH_TO_SPEX_UTIL) -O
I = -I../Include -I../Source -I../../../SuiteSparse_config -I../../../COLAMD/Include -I../../../AMD/Include -I../../SPEX_Util/Include -I../../SPEX_Util/Source -I../../SPEX_Util/Lib

LDLIBS += -lm -lsuitesparseconfig -lgmp -lmpfr -lcolamd -lamd -L../../SPEX_Util/Lib/ -lspexutil 
CS = ../Lib/libspexwamf.a $(LDLIBS)


all: lib example
	- ./example

lib:
	( cd ../Lib ; $(MAKE) )

example: lib example.c Makefile
	$(CC) $(LDFLAGS) $(CF) $(I) -g -o example example.c demos.c $(CS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r SPEX_Chol SPEX_Chol_debug *.a *.dSYM *.obj *.dll
