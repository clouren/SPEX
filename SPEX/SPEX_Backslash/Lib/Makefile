#-------------------------------------------------------------------------------
# SPEX_Backslash/Lib/Makefile
#-------------------------------------------------------------------------------

# To just compile the library, do:
#       make
#
# To run a demo using the library
#       cd ../Demo ; make
#


LIBRARY = libspexbackslash
VERSION = 1.0.1
SO_VERSION = 1

default: library

SUITESPARSE ?= $(realpath $(CURDIR)/../../..)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

PATH_TO_SPEX_UTIL=-Wl,-rpath,"../../SPEX_Util/Lib/"
PATH_TO_SPEX_CHOL=-Wl,-rpath,"../../SPEX_Cholesky/Lib/"
PATH_TO_SPEX_LEFT=-Wl,-rpath,"../../SPEX_Left_LU/Lib/"
PATH_TO_SPEX_BACK=-Wl,-rpath,"../SPEX_Backslash/Lib/"


# uncomment for extra error checking:
#CFLAGS += -Wall -Wextra -Wpedantic -Werror

# SPEX_Backslash depends on SuiteSparse_config, AMD, COLAMD, M, GMP, MPFR, SPEX_Util
# SPEX_Left_LU and SPEX_Cholesky
LDLIBS += -lsuitesparseconfig -lamd -lcolamd \
-L../../SPEX_Util/Lib/ -lspexutil \
-L../../SPEX_Left_LU/Lib/ -lspexleftlu \
-L../../SPEX_Cholesky/Lib/ -lspexchol \
-lm -lgmp -lmpfr
CF += -g
C = $(CC) $(CF) -I../Include -I../../../COLAMD/Include -I../../../AMD/Include -I../../../SuiteSparse_config -I../../../include/ -I../../SPEX_Util/Include -I../../SPEX_Util/Source -I../../SPEX_Left_LU/Include -I../../SPEX_Left_LU/Source -I../../SPEX_Cholesky/Include -I../../SPEX_Cholesky/Source
#C += -g

all: install

OBJ = \
     SPEX_Backslash.o
$(OBJ): ../Include/SPEX_Backslash.h Makefile

%.o: ../Source/%.c ../Include/SPEX_Backslash.h 
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
install: $(AR_TARGET) $(INSTALL_LIB)/$(SO_TARGET)

$(INSTALL_LIB)/$(SO_TARGET): $(OBJ)
	@mkdir -p $(INSTALL_LIB)
	@mkdir -p $(INSTALL_INCLUDE)
	@mkdir -p $(INSTALL_DOC)
	$(CC) $(SO_OPTS) $^ -o $@ $(LDLIBS)
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_PLAIN) )
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_MAIN) )
	$(CP) ../Include/SPEX_Backslash.h $(INSTALL_INCLUDE)
#	$(CP) ../Doc/SPEX_LU_UserGuide.pdf $(INSTALL_DOC)
#	$(CP) ../README.md $(INSTALL_DOC)/SPEX_LU_README.md
	chmod 755 $(INSTALL_LIB)/$(SO_TARGET)
	chmod 644 $(INSTALL_INCLUDE)/SPEX_Backslash.h
#	chmod 644 $(INSTALL_DOC)/SPEX_LU_UserGuide.pdf
#	chmod 644 $(INSTALL_DOC)/SPEX_LU_README.md

# uninstall SPEX_Chol
uninstall:
	$(RM) $(INSTALL_LIB)/$(SO_TARGET)
	$(RM) $(INSTALL_LIB)/$(SO_PLAIN)
	$(RM) $(INSTALL_LIB)/$(SO_MAIN)
	$(RM) $(INSTALL_INCLUDE)/SPEX_Backslash.h
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
