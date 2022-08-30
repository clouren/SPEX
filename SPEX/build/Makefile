#-------------------------------------------------------------------------------
# SPEX/build/Makefile
#
# SPEX: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
# Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.
# All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later
#-------------------------------------------------------------------------------

# To just compile the library, do:
#       make
#
# To run a demo using the library
#       cd ../Demo ; make

LIBRARY = libspex
VERSION = 1.1.0
SO_VERSION = 1

default: library

SUITESPARSE ?= $(realpath $(CURDIR)/../..)
include ../../SuiteSparse_config/SuiteSparse_config.mk

# uncomment for extra error checking:
CFLAGS += -Wall -Wextra -Wpedantic -Werror

# SPEX depends on SuiteSparse_config, AMD, COLAMD, M, GMP, and MPFR
LDLIBS += -lsuitesparseconfig -lamd -lcolamd -lm -lgmp -lmpfr

# uncomment for debugging
#CF += -g

C = $(CC) $(CF) -I../SPEX_Backslash/Include -I../../COLAMD/Include -I../../AMD/Include -I../../SuiteSparse_config -I../SPEX_Utilities/Source -I../SPEX_LU/Source -I../SPEX_Cholesky/Source -I../SPEX_Update/Source -I../SPEX_Backslash/Source

all: install

#-------------------------------------------------------------------------------

OBJ_UTIL = \
    SPEX_calloc.o \
    spex_cast_array.o \
    spex_cast_matrix.o \
    SPEX_create_default_options.o \
    spex_create_mpfr_array.o \
    spex_create_mpq_array.o \
    spex_create_mpz_array.o \
    spex_cumsum.o \
    spex_expand_double_array.o \
    spex_expand_mpfr_array.o \
    spex_expand_mpq_array.o \
    SPEX_finalize.o \
    SPEX_free.o \
    SPEX_gmp.o \
    SPEX_initialize.o \
    SPEX_initialize_expert.o \
    SPEX_malloc.o \
    SPEX_matrix_allocate.o \
    SPEX_matrix_check.o \
    SPEX_matrix_copy.o \
    SPEX_matrix_free.o \
    spex_matrix_mul.o \
    SPEX_matrix_nnz.o \
    spex_permute_dense_matrix.o \
    SPEX_realloc.o \
    spex_sparse_collapse.o \
    spex_sparse_realloc.o \
    SPEX_symbolic_analysis_free.o \
    SPEX_factorization_free.o \
    SPEX_vector_allocate.o \
    SPEX_vector_realloc.o \
    SPEX_vector_free.o \
    spex_dynamic_to_CSC_mpz.o \
    spex_CSC_mpz_to_dynamic.o \
    spex_delete_mpz_array.o \
    SPEX_determine_symmetry.o SPEX_transpose.o \
    SPEX_factorization_convert.o \
    SPEX_factorization_check.o \
    spex_factorization_basic_check.o \
    spex_colamd.o \
    spex_amd.o
    # spex_matrix_div.o \

$(OBJ_UTIL): ../SPEX_Backslash/Include/SPEX.h ../SPEX_Utilities/Source/spex_util_internal.h Makefile

SPEX_%.o: ../SPEX_Utilities/Source/SPEX_%.c ../SPEX_Backslash/Include/SPEX.h ../SPEX_Utilities/Source/spex_util_internal.h
	$(C) $(I) -c $<

spex_%.o: ../SPEX_Utilities/Source/spex_%.c ../SPEX_Backslash/Include/SPEX.h ../SPEX_Utilities/Source/spex_util_internal.h
	$(C) $(I) -c $<

#-------------------------------------------------------------------------------

OBJ_LEFTLU = \
    SPEX_lu_analyze.o \
    SPEX_lu_backslash.o \
    spex_left_lu_back_sub.o \
    spex_left_lu_dfs.o \
    SPEX_lu_factorize.o \
    spex_left_lu_forward_sub.o \
    spex_left_lu_get_largest_pivot.o \
    spex_left_lu_get_nonzero_pivot.o \
    spex_left_lu_get_pivot.o \
    spex_left_lu_get_smallest_pivot.o \
    spex_left_lu_reach.o \
    spex_left_lu_ref_triangular_solve.o \
    SPEX_lu_solve.o

$(OBJ_LEFTLU): ../SPEX_Backslash/Include/SPEX.h ../SPEX_LU/Source/spex_lu_internal.h Makefile

SPEX_lu_%.o: ../SPEX_LU/Source/SPEX_lu_%.c ../SPEX_Backslash/Include/SPEX.h ../SPEX_LU/Source/spex_lu_internal.h
	$(C) $(I) -c $<

SPEX_lu_%.o: ../SPEX_LU/Source/SPEX_lu_%.c ../SPEX_Backslash/Include/SPEX.h ../SPEX_LU/Source/spex_lu_internal.h
	$(C) $(I) -c $<

spex_left_lu_%.o: ../SPEX_LU/Source/spex_left_lu_%.c ../SPEX_Backslash/Include/SPEX.h ../SPEX_LU/Source/spex_lu_internal.h
	$(C) $(I) -c $<

#-------------------------------------------------------------------------------

OBJ_CHOL = \
    SPEX_cholesky_analyze.o \
    SPEX_cholesky_backslash.o \
    spex_cholesky_backward_sub.o \
    spex_cholesky_counts.o \
    spex_cholesky_ereach.o \
    spex_cholesky_etree.o \
    spex_cholesky_factor.o \
    SPEX_cholesky_factorize.o \
    spex_cholesky_forward_sub.o \
    spex_cholesky_leaf.o \
    spex_cholesky_left_factor.o \
    spex_cholesky_permute_A.o \
    spex_cholesky_post.o \
    spex_cholesky_pre_left_factor.o \
    spex_cholesky_preorder.o \
    SPEX_cholesky_solve.o \
    spex_cholesky_symbolic_analysis.o \
    spex_cholesky_tdfs.o \
    spex_cholesky_up_factor.o \
    spex_cholesky_left_triangular_solve.o \
    spex_cholesky_up_triangular_solve.o

$(OBJ_CHOL): ../SPEX_Backslash/Include/SPEX.h ../SPEX_Cholesky//Source/spex_cholesky_internal.h ../SPEX_Utilities/Source/spex_util_internal.h Makefile

SPEX_cholesky_%.o: ../SPEX_Cholesky/Source/SPEX_cholesky_%.c ../SPEX_Backslash/Include/SPEX.h ../SPEX_Cholesky/Source/spex_cholesky_internal.h
	$(C) $(I) -c $<

spex_cholesky_%.o: ../SPEX_Cholesky/Source/spex_cholesky_%.c ../SPEX_Backslash/Include/SPEX.h ../SPEX_Cholesky/Source/spex_cholesky_internal.h
	$(C) $(I) -c $<

#-------------------------------------------------------------------------------

OBJ_UPDATE = \
    spex_update_backward_sub.o \
    SPEX_update_cholesky_rank1.o \
    spex_update_cppu.o \
    spex_update_dppu1.o \
    spex_update_dppu2.o \
    spex_update_finalize_and_insert_vk.o \
    spex_update_forward_sub.o \
    spex_update_get_scattered_v.o \
    spex_update_insert_new_entry.o \
    spex_update_ipge.o \
    SPEX_update_lu_colrep.o \
    SPEX_update_matrix_colrep.o \
    SPEX_update_solve.o \
    SPEX_update_tsolve.o \
    spex_update_solve_internal.o \
    spex_update_triangular_solve.o \
    spex_update_debug.o

$(OBJ_UPDATE): ../SPEX_Backslash/Include/SPEX.h ../SPEX_Update//Source/spex_update_internal.h Makefile

SPEX_update_%.o: ../SPEX_Update/Source/SPEX_update_%.c ../SPEX_Backslash/Include/SPEX.h ../SPEX_Update/Source/spex_update_internal.h
	$(C) $(I) -c $<

spex_update_%.o: ../SPEX_Update/Source/spex_update_%.c ../SPEX_Backslash/Include/SPEX.h ../SPEX_Update/Source/spex_update_internal.h
	$(C) $(I) -c $<

#-------------------------------------------------------------------------------

OBJ_BACKSLASH = \
    SPEX_backslash.o

SPEX_backslash.o: ../SPEX_Backslash/Source/SPEX_backslash.c ../SPEX_Backslash/Include/SPEX.h
	$(C) $(I) -c $<

#-------------------------------------------------------------------------------
# Create the bin
#-------------------------------------------------------------------------------

library:
	$(MAKE) install INSTALL=$(SUITESPARSE)

static: $(AR_TARGET)

$(AR_TARGET): $(OBJ_UTIL) $(OBJ_LEFTLU) $(OBJ_UPDATE) $(OBJ_CHOL) $(OBJ_BACKSLASH)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $(AR_TARGET)

#-------------------------------------------------------------------------------
# install (shared C version only)
#-------------------------------------------------------------------------------

# install SPEX
install: $(AR_TARGET) $(INSTALL_LIB)/$(SO_TARGET)

$(INSTALL_LIB)/$(SO_TARGET): $(OBJ_UTIL) $(OBJ_LEFTLU) $(OBJ_UPDATE) $(OBJ_CHOL) $(OBJ_BACKSLASH)
	@mkdir -p $(INSTALL_LIB)
	@mkdir -p $(INSTALL_INCLUDE)
	@mkdir -p $(INSTALL_DOC)
	$(CC) $(SO_OPTS) $^ -o $@ $(LDLIBS)
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_PLAIN) )
	( cd $(INSTALL_LIB) ; ln -sf $(SO_TARGET) $(SO_MAIN) )
	$(CP) ../SPEX_Backslash/Include/SPEX.h $(INSTALL_INCLUDE)
	#$(CP) ../Doc/SPEX_UserGuide.pdf $(INSTALL_DOC)
	#$(CP) ../README.md $(INSTALL_DOC)/SPEX_README.md
	chmod 644 $(INSTALL_INCLUDE)/SPEX.h
	#chmod 644 $(INSTALL_DOC)/SPEX_UserGuide.pdf
	#chmod 644 $(INSTALL_DOC)/SPEX_README.md


# uninstall SPEX
uninstall:
	$(RM) $(INSTALL_LIB)/$(SO_TARGET)
	$(RM) $(INSTALL_LIB)/$(SO_PLAIN)
	$(RM) $(INSTALL_LIB)/$(SO_MAIN)
	$(RM) $(INSTALL_INCLUDE)/SPEX.h
	$(RM) $(INSTALL_DOC)/SPEX_UserGuide.pdf
	$(RM) $(INSTALL_DOC)/SPEX_README.md


#-------------------------------------------------------------------------------
# Remove all but the files in the original distribution
#-------------------------------------------------------------------------------

clean:
	- $(RM) -r $(CLEAN)

purge: distclean

distclean: clean
	- $(RM) -r $(PURGE)

