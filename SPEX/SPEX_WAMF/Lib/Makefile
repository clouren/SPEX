# On the Mac, you need to first install Xcode.  Then
# type this command in the Terminal:
#       xcode-select --install
# to ensure that you have the command-line tools installed.
#
# To just compile the library, do:
#       make
#
# To run a demo using the library 
#       cd ../Demo ; make
#

LIBRARY = libspexwamf
VERSION = 1.0.1
SO_VERSION = 1

default: library

SUITESPARSE ?= $(realpath $(CURDIR)/../../..)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
#CFLAGS += -Wall -Wextra -Wpedantic -Werror

# SPEX_Chol depends on SuiteSparse_config, AMD, COLAMD, M, GMP, and MPFR
LDLIBS += -lsuitesparseconfig -lamd -lcolamd -lm -lgmp -lmpfr

C = $(CC) $(CF) -I../Include -I../../../COLAMD/Include -I../../../AMD/Include -I../../../SuiteSparse_config -I../../../include/ -I../../SPEX_Util/Include -I../../SPEX_Util/Source

all: install

OBJ = \
     SPEX_WAMF_diag.o SPEX_WAMF_scatter.o SPEX_WAMF_tdfs.o SPEX_WAMF_fkeep.o \
     SPEX_WAMF_sparse_realloc.o SPEX_WAMF_wamf.o SPEX_WAMF_get_bitmat.o \
     SPEX_WAMF_symbolic_add.o SPEX_WAMF_wclear.o SPEX_WAMF_get_weights.o \
     SPEX_WAMF_symbolic_multiply.o
     

$(OBJ): ../Include/SPEX_WAMF.h Makefile

%.o: ../Source/%.c ../Include/SPEX_WAMF.h 
	$(C) $(I) -c $<

#-------------------------------------------------------------------------------
# Create the bin
#-------------------------------------------------------------------------------

library:
	$(MAKE) install INSTALL=$(SUITESPARSE)

static: $(AR_TARGET)

$(AR_TARGET): $(OBJ)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $(AR_TARGET)

#-------------------------------------------------------------------------------
# install (shared C version only)
#-------------------------------------------------------------------------------

# install SPEX_LU
install: $(AR_TARGET) #$(INSTALL_LIB)/$(SO_TARGET)

$(INSTALL_LIB)/$(SO_TARGET): $(OBJ)
	@mkdir -p $(INSTALL_LIB)
	@mkdir -p $(INSTALL_INCLUDE)
	@mkdir -p $(INSTALL_DOC)
	$(CC) $(SO_OPTS) $^ -o $@ $(LDLIBS)
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_PLAIN) )
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_MAIN) )
	$(CP) ../Include/SPEX_WAMF.h $(INSTALL_INCLUDE)
#	$(CP) ../Doc/SPEX_LU_UserGuide.pdf $(INSTALL_DOC)
#	$(CP) ../README.md $(INSTALL_DOC)/SPEX_LU_README.md
	chmod 755 $(INSTALL_LIB)/$(SO_TARGET)
	chmod 644 $(INSTALL_INCLUDE)/SPEX_WAMF.h
#	chmod 644 $(INSTALL_DOC)/SPEX_LU_UserGuide.pdf
#	chmod 644 $(INSTALL_DOC)/SPEX_LU_README.md

# uninstall SPEX_Chol
uninstall:
	$(RM) $(INSTALL_LIB)/$(SO_TARGET)
	$(RM) $(INSTALL_LIB)/$(SO_PLAIN)
	$(RM) $(INSTALL_LIB)/$(SO_MAIN)
	$(RM) $(INSTALL_INCLUDE)/SPEX_WAMF.h
#	$(RM) $(INSTALL_DOC)/SPEX_LU_UserGuide.pdf
#	$(RM) $(INSTALL_DOC)/SPEX_LU_README.md

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) -r $(CLEAN)

purge: distclean

distclean: clean
	- $(RM) -r $(PURGE)
