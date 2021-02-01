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

include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
#CFLAGS += -Wall -Wextra -Wpedantic -Werror

# SPEX_Chol depends on SuiteSparse_config, AMD, COLAMD, M, GMP, and MPFR
LDLIBS += -lsuitesparseconfig -lamd -lcolamd -lm -lgmp -lmpfr

C = $(CC) $(CF) -I../Include -I../../COLAMD/Include -I../../AMD/Include -I../../SuiteSparse_config -I../../include/ -I../../SLIP_Util/Include -I../../SLIP_Util/Source

all: install

OBJ = \
     SPEX_Chol_analysis_free.o spex_Chol_post.o spex_Chol_counts.o spex_Chol_Pre_Left_Factor.o \
     spex_Chol_ereach.o \
     SPEX_Chol_Solve.o spex_Chol_etree.o spex_Chol_tdfs.o SPEX_Chol_Factor.o \
     spex_Chol_forward_sub.o spex_Left_Chol_triangular_solve.o \
     spex_Chol_leaf.o spex_Chol_ltsolve.o SPEX_Chol_permute_A.o \
     spex_Up_Chol_triangular_solve.o  SPEX_Chol_backslash.o \
     SPEX_Chol_preorder.o 

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
install: $(AR_TARGET) #$(INSTALL_LIB)/$(SO_TARGET)

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
