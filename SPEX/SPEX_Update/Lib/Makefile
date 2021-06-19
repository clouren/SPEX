#-------------------------------------------------------------------------------
# SPEX_Update/Lib/Makefile
#-------------------------------------------------------------------------------

# To just compile the library, do:
#       make
#
# To run a demo using the library
#       cd ../Demo ; make

LIBRARY = libspexupdate
VERSION = 1.0.0
SO_VERSION = 1

default: library

SUITESPARSE ?= $(realpath $(CURDIR)/../../..)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
CFLAGS += -Wall -Wextra -Wpedantic -Werror

# SPEX_CHOLMOD depends on SuiteSparse_config, M, GMP, and MPFR
LDLIBS += -lspexutil -lsuitesparseconfig -lm -lgmp -lmpfr

C = $(CC) $(CF) -I../Include -I../../../SuiteSparse_config -I../../../include/ -I../../SPEX_Util/Include -I../../SPEX_Util/Source

all: install

OBJ = \
    spex_update_backward_sub.o             spex_update_get_scattered_v.o \
    spex_update_cppu.o \
    spex_update_insert_new_entry.o \
    spex_update_ipge.o \
    SPEX_Update_LU_ColRep.o \
    spex_update_dppu1.o\
    spex_update_dppu2.o\
    spex_update_finalize_and_insert_vk.o   SPEX_Update_Solve.o \
    spex_update_triangular_solve.o \
    spex_update_forward_sub.o \
    spex_update_verify.o \
    SPEX_Update_Chol_Rank1.o \
    SPEX_Update_permute_row.o SPEX_Update_matrix_canonicalize.o


$(OBJ): ../Include/SPEX_Update.h ../Source/spex_update_internal.h Makefile

%.o: ../Source/%.c ../Include/SPEX_Update.h ../Source/spex_update_internal.h
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

# install SPEX_Update
install: $(AR_TARGET) $(INSTALL_LIB)/$(SO_TARGET)

$(INSTALL_LIB)/$(SO_TARGET): $(OBJ)
	@mkdir -p $(INSTALL_LIB)
	@mkdir -p $(INSTALL_INCLUDE)
	@mkdir -p $(INSTALL_DOC)
	$(CC) $(SO_OPTS) $^ -o $@ $(LDLIBS)
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_PLAIN) )
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_MAIN) )
	$(CP) ../Include/SPEX_Update.h $(INSTALL_INCLUDE)
	#$(CP) ../Doc/SLIP_LU_UserGuide.pdf $(INSTALL_DOC)
	#$(CP) ../README.md $(INSTALL_DOC)/SLIP_LU_README.md
	chmod 755 $(INSTALL_LIB)/$(SO_TARGET)
	chmod 644 $(INSTALL_INCLUDE)/SPEX_Update.h
	#chmod 644 $(INSTALL_DOC)/SLIP_LU_UserGuide.pdf
	#chmod 644 $(INSTALL_DOC)/SLIP_LU_README.md

# uninstall SPEX_Update
uninstall:
	$(RM) $(INSTALL_LIB)/$(SO_TARGET)
	$(RM) $(INSTALL_LIB)/$(SO_PLAIN)
	$(RM) $(INSTALL_LIB)/$(SO_MAIN)
	$(RM) $(INSTALL_INCLUDE)/SPEX_Update.h
	$(RM) $(INSTALL_DOC)/SLIP_LU_UserGuide.pdf
	$(RM) $(INSTALL_DOC)/SLIP_LU_README.md

#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) -r $(CLEAN)

purge: distclean

distclean: clean
	- $(RM) -r $(PURGE)

