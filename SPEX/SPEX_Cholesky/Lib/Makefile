#-------------------------------------------------------------------------------
# SPEX_Chol/Lib/Makefile
#-------------------------------------------------------------------------------

# To just compile the library, do:
#       make
#
# To run a demo using the library
#       cd ../Demo ; make
#


LIBRARY = libspexchol
VERSION = 1.0.1
SO_VERSION = 1

default: library

SUITESPARSE ?= $(realpath $(CURDIR)/../../..)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
#CFLAGS += -Wall -Wextra -Wpedantic -Werror

# SPEX_Chol depends on SuiteSparse_config, AMD, COLAMD, M, GMP, and MPFR
LDLIBS += -lsuitesparseconfig -lamd -lcolamd -lspexutil -lm -lgmp -lmpfr
CF += -g
C = $(CC) $(CF) -I../Include -I../../../COLAMD/Include -I../../../AMD/Include -I../../../SuiteSparse_config -I../../../include/ -I../../SPEX_Util/Include -I../../SPEX_Util/Source
#C += -g

all: install

OBJ = \
     spex_chol_post.o spex_chol_counts.o spex_chol_pre_left_factor.o \
     spex_chol_ereach.o \
     SPEX_Chol_solve.o spex_chol_etree.o spex_chol_tdfs.o spex_chol_factor.o \
     spex_chol_forward_sub.o spex_left_chol_triangular_solve.o \
     spex_chol_leaf.o spex_chol_backward_sub.o spex_chol_permute_A.o \
     spex_up_chol_triangular_solve.o  SPEX_Chol_backslash.o \
     spex_chol_preorder.o spex_chol_up_factor.o spex_chol_left_factor.o \
     spex_chol_symbolic_analysis.o SPEX_Chol_analyze.o\
     SPEX_Chol_factorize.o \

$(OBJ): ../Include/SPEX_Chol.h Makefile

%.o: ../Source/%.c ../Include/SPEX_Chol.h 
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
	$(CP) ../Include/SPEX_Chol.h $(INSTALL_INCLUDE)
#	$(CP) ../Doc/SPEX_LU_UserGuide.pdf $(INSTALL_DOC)
#	$(CP) ../README.md $(INSTALL_DOC)/SPEX_LU_README.md
	chmod 755 $(INSTALL_LIB)/$(SO_TARGET)
	chmod 644 $(INSTALL_INCLUDE)/SPEX_Chol.h
#	chmod 644 $(INSTALL_DOC)/SPEX_LU_UserGuide.pdf
#	chmod 644 $(INSTALL_DOC)/SPEX_LU_README.md

# uninstall SPEX_Chol
uninstall:
	$(RM) $(INSTALL_LIB)/$(SO_TARGET)
	$(RM) $(INSTALL_LIB)/$(SO_PLAIN)
	$(RM) $(INSTALL_LIB)/$(SO_MAIN)
	$(RM) $(INSTALL_INCLUDE)/SPEX_Chol.h
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
