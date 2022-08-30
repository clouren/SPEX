#-------------------------------------------------------------------------------
# SPEX/SPEX/Tcov/Makefile: compile and run SLIP test coverage
#-------------------------------------------------------------------------------

# SPEX: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
# Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#-------------------------------------------------------------------------------

default: cov

SUITESPARSE ?= $(realpath $(CURDIR)/../..)
include ../../SuiteSparse_config/SuiteSparse_config.mk

# Linux test coverage (gcc is required for test coverage)
CC = gcc
CFLAGS = -g -fprofile-arcs -ftest-coverage \
        -Wall -W -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
        -Wredundant-decls -Wnested-externs -Wdisabled-optimization -std=c99 \
        -Wno-unused-parameter \
	-I../SPEX_Backslash/Include/ \
	-I../SPEX_Utilities/Source/ \
	-I../SPEX_Cholesky/Source/ -I../SPEX_Cholesky/Demo/ \
	-I../SPEX_LU/Source -I../SPEX_LU/Demo \
	-I../SPEX_Update/Source \
        -I../../SuiteSparse_config -I../../COLAMD/Include \
	-I../../AMD/Include \
        -DSPEX_GMP_LIST_INIT=2 #-DSPEX_DEBUG

LDLIBS += -lm -lgmp -lmpfr -lcolamd -lamd -lsuitesparseconfig

# run all statement coverage tests, and then check for 100% coverage
cov: runtests
	./covall

all: test4llu SPEXLLU test4luu test4chol

#-------------------------------------------------------------------------------
# compile .c file in this folder

OBJ_Tcov = tcov_malloc_test.o

$(OBJ_Tcov): tcov_malloc_test.h

#-------------------------------------------------------------------------------
# link and compile each .c file in SPEX_Utilities/Source

OBJ_Util = \
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
    SPEX_factorization_convert.o \
    SPEX_factorization_check.o \
    spex_factorization_basic_check.o \
    SPEX_determine_symmetry.o \
    SPEX_transpose.o \
    spex_amd.o \
    spex_colamd.o

$(OBJ_Util): ../SPEX_Backslash/Include/SPEX.h ../SPEX_Utilities/Source/spex_util_internal.h

.PRECIOUS: SPEX_%.c spex_%.c

SPEX_%.c:
	- ln -s ../SPEX_Utilities/Source/$@

spex_%.c:
	- ln -s ../SPEX_Utilities/Source/$@

#-------------------------------------------------------------------------------
# link and compile each .c file in SPEX_LU/Source

OBJ_LU = \
    SPEX_LU_analyze.o \
    SPEX_LU_backslash.o \
    spex_left_lu_back_sub.o \
    spex_left_lu_dfs.o \
    SPEX_LU_factorize.o \
    spex_left_lu_forward_sub.o \
    spex_left_lu_get_largest_pivot.o \
    spex_left_lu_get_nonzero_pivot.o \
    spex_left_lu_get_pivot.o \
    spex_left_lu_get_smallest_pivot.o \
    spex_left_lu_reach.o \
    spex_left_lu_ref_triangular_solve.o \
    SPEX_LU_solve.o

.PRECIOUS: SPEX_LU_%.c SPEX_LU_%.c spex_left_lu_%.c spexlu_demo.c demos.c

SPEX_LU_%.c:
	- ln -s ../SPEX_LU/Source/$@

SPEX_LU_%.c:
	- ln -s ../SPEX_LU/Source/$@

spex_left_lu_%.c:
	- ln -s ../SPEX_LU/Source/$@

spexlu_demo.c:
	- ln -s ../SPEX_LU/Demo/$@

demos.c:
	- ln -s ../SPEX_LU/Demo/$@

#-------------------------------------------------------------------------------
# link and compile each .c file in SPEX_Update/Source

OBJ_Update = \
    spex_update_backward_sub.o \
    SPEX_Update_Chol_Rank1.o \
    spex_update_cppu.o \
    spex_update_dppu1.o \
    spex_update_dppu2.o \
    spex_update_finalize_and_insert_vk.o \
    spex_update_forward_sub.o \
    spex_update_get_scattered_v.o \
    spex_update_insert_new_entry.o \
    spex_update_ipge.o \
    SPEX_Update_LU_ColRep.o \
    SPEX_Update_matrix_colrep.o \
    SPEX_Update_solve.o \
    SPEX_Update_tsolve.o \
    spex_update_solve_internal.o \
    spex_update_triangular_solve.o


$(OBJ_Update): ../SPEX_Backslash/Include/SPEX.h ../SPEX_Update/Source/spex_update_internal.h

.PRECIOUS: SPEX_Update_%.c spex_update_%.c

SPEX_Update_%.c:
	- ln -s ../SPEX_Update/Source/$@

spex_update_%.c:
	- ln -s ../SPEX_Update/Source/$@

#-------------------------------------------------------------------------------
# link and compile each .c file in SPEX_Cholesky/Source

OBJ_Cholesky = \
    SPEX_Chol_analyze.o \
    SPEX_Chol_backslash.o \
    spex_chol_backward_sub.o \
    spex_chol_counts.o \
    spex_chol_ereach.o \
    spex_chol_etree.o \
    spex_chol_factor.o \
    SPEX_Chol_factorize.o \
    spex_chol_forward_sub.o \
    spex_chol_leaf.o \
    spex_chol_left_factor.o \
    spex_chol_left_triangular_solve.o \
    spex_chol_permute_A.o \
    spex_chol_post.o \
    spex_chol_pre_left_factor.o \
    spex_chol_preorder.o \
    SPEX_Chol_solve.o \
    spex_chol_symbolic_analysis.o \
    spex_chol_tdfs.o \
    spex_chol_up_factor.o \
    spex_chol_up_triangular_solve.o

$(OBJ_Cholesky): ../SPEX_Backslash/Include/SPEX.h

.PRECIOUS: SPEX_Chol_%.c spex_chol_%.c

SPEX_Chol_%.c:
	- ln -s ../SPEX_Cholesky/Source/$@

spex_chol_%.c:
	- ln -s ../SPEX_Cholesky/Source/$@

chol_demos.c:
	- ln -s ../SPEX_Cholesky/Demo/$@

#-------------------------------------------------------------------------------

test4llu: $(OBJ_Tcov) $(OBJ_Util) $(OBJ_LU) $(OBJ_Update) tcov_for_left_LU.c
	$(CC) $(LDFLAGS) tcov_for_left_LU.c $(CFLAGS) -o tcov_for_left_LU $(OBJ_Tcov) $(OBJ_Util) $(OBJ_LU) $(OBJ_Update) $(LDLIBS)

SPEXLLU: $(OBJ_Util) $(OBJ_LU) $(OBJ_Update) spexlu_demo.c demos.c ../SPEX_LU/Demo/demos.h
	$(CC) $(LDFLAGS) spexlu_demo.c demos.c $(CFLAGS) -o spexlu_demo $(OBJ_Tcov) $(OBJ_Util) $(OBJ_LU) $(OBJ_Update) $(LDLIBS)

test4luu: $(OBJ_Tcov) $(OBJ_Util) $(OBJ_LU) $(OBJ_Update) tcov_for_update.c simple_rand.c simple_rand.h
	$(CC) $(LDFLAGS) tcov_for_update.c simple_rand.c $(CFLAGS) -o tcov_for_update $(OBJ_Tcov) $(OBJ_Util) $(OBJ_LU) $(OBJ_Update) $(LDLIBS)

test4chol: $(OBJ_Tcov) $(OBJ_Util) $(OBJ_Cholesky) tcov_for_cholesky.c simple_rand.c simple_rand.h chol_demos.c
	$(CC) $(LDFLAGS) $(CFLAGS) tcov_for_cholesky.c simple_rand.c chol_demos.c -o tcov_for_cholesky $(OBJ_Tcov) $(OBJ_Util) $(OBJ_Cholesky) $(LDLIBS)

#-------------------------------------------------------------------------------

# run all statement coverage tests
runtests: runtests4llu runtests4luu runtests4chol

# only run test for SPEX_LU
runtests4llu: all
	- ./spexlu_demo p 2 q 1 f ../SPEX_LU/ExampleMats/10teams_mat.txt ../SPEX_LU/ExampleMats/10teams_v.txt
	- ./spexlu_demo p 3 q 2 o 1 f ../SPEX_LU/ExampleMats/10teams_mat.txt ../SPEX_LU/ExampleMats/10teams_v.txt
	- ./spexlu_demo p 4 q 3 o 1 f ../SPEX_LU/ExampleMats/test_mat.txt ../SPEX_LU/ExampleMats/test_rhs.txt
	- ./spexlu_demo p 5 f ../SPEX_LU/ExampleMats/10teams_mat.txt ../SPEX_LU/ExampleMats/10teams_v.txt
	- ./tcov_for_left_LU > test_for_left_LU.out
	- ./tcov_for_left_LU 0 1 1

# only run test for SPEX_Update
runtests4luu: all
	- ./tcov_for_update 0 6 > test_for_update1.out
	- ./tcov_for_update 0 0 26 > test_for_update2.out
	- ./tcov_for_update 1 6 > test_for_update3.out

# only run test for SPEX_Cholesky
runtests4chol: test4chol
	- ./tcov_for_cholesky ../SPEX_Cholesky/ExampleMats/872.mat.txt

#-------------------------------------------------------------------------------

# To run with valgrind:
V = valgrind #--leak-check=full

# run all statement coverage tests but with valgrind
vtests: vtests4llu vtests4luu vtests4chol

# run test for SPEX_LU with valgrind
vtests4llu: all
	- $(V) ./spexlu_demo p 2 q 0 f ../SPEX_LU/ExampleMats/10teams_mat.txt ../SPEX_LU/ExampleMats/10teams_v.txt
	- $(V) ./spexlu_demo p 3 q 1 o 1 f ../SPEX_LU/ExampleMats/10teams_mat.txt ../SPEX_LU/ExampleMats/10teams_v.txt
	- $(V) ./spexlu_demo p 4 q 2 o 1 f ../SPEX_LU/ExampleMats/test_mat.txt ../SPEX_LU/ExampleMats/test_rhs.txt
	- $(V) ./spexlu_demo p 5 f ../SPEX_LU/ExampleMats/10teams_mat.txt ../SPEX_LU/ExampleMats/10teams_v.txt
	- $(V) ./tcov_for_left_LU > test_for_left_LU.out
	- $(V) ./tcov_for_left_LU 0 1 1

# run test for SPEX_Update with valgrind
vtests4luu: all
	- $(V) ./tcov_for_update 0 6 > test_for_update1.out
	- $(V) ./tcov_for_update 0 0 25 > test_for_update2.out
	- $(V) ./tcov_for_update 1 6 > test_for_update3.out

# run test for SPEX_Cholesky with valgrind
vtests4chol: test4chol
	- $(V) --leak-check=full ./tcov_for_cholesky ../SPEX_Cholesky/ExampleMats/872.mat.txt

#-------------------------------------------------------------------------------

# remove all files not in the original distribution, including symbolic links
clean:
	- $(RM) *.o *.bbg *.da *.gcov *.gcda *gcno
	- $(RM) tcov_for_left_LU tcov_for_update spexlu_demo *.out *.a cov.sort out Mats4Tcov/SPEX_Update/mat.txt
	- $(RM) -r SPEX_*.c spex_*.c *.dSYM spexlu_demo.c demos.c
	- $(RM) -r tcov_for_cholesky chol_demos.c

purge: distclean

distclean: clean
