#-------------------------------------------------------------------------------
# SPEX/SPEX/Tcov/Makefile: compile and run SLIP test coverage
#-------------------------------------------------------------------------------

# SPEX: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
# Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#-------------------------------------------------------------------------------

default: run

SUITESPARSE ?= $(realpath $(CURDIR)/../..)
include ../../SuiteSparse_config/SuiteSparse_config.mk

# Linux test coverage (gcc is required for test coverage)
CC = gcc
CFLAGS = -g -fprofile-arcs -ftest-coverage \
        -Wall -W -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
        -Wredundant-decls -Wnested-externs -Wdisabled-optimization -std=c99 \
        -Wno-unused-parameter \
	-I../SPEX_Util/Include/ -I../SPEX_Util/Source/ \
	-I../SPEX_Left_LU/Include -I../SPEX_Left_LU/Source \
	-I../SPEX_Left_LU/Demo \
	-I../SPEX_Update/Include -I../SPEX_Update/Source \
	-I../SPEX_Update/Demo \
        -I../../SuiteSparse_config -I../../COLAMD/Include \
	-I../../AMD/Include \
        -DSPEX_GMP_LIST_INIT=2 #-DSPEX_DEBUG

LDLIBS += -lm -lgmp -lmpfr -lcolamd -lamd -lsuitesparseconfig

# run all statement coverage tests, and then check for 100% coverage
run: runtests
	./covall

all: test4llu SPEXLLU test4luu SPEXLUU

CS_Tcov = tcov_malloc_test.o

$(CS_Tcov): tcov_malloc_test.h

CS_Util = \
    SPEX_calloc.o \
    spex_cast_array.o \
    spex_cast_matrix.o \
    SPEX_check_solution.o \
    SPEX_create_default_options.o \
    spex_create_mpfr_array.o \
    spex_create_mpq_array.o \
    spex_create_mpz_array.o \
    SPEX_cumsum.o \
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
    SPEX_matrix_div.o \
    SPEX_matrix_free.o \
    SPEX_matrix_mul.o \
    SPEX_matrix_nnz.o \
    SPEX_realloc.o \
    spex_sparse_collapse.o \
    spex_sparse_realloc.o \
    SPEX_LU_analysis_free.o \
    SPEX_vector_allocate.o \
    SPEX_vector_realloc.o \
    SPEX_vector_free.o \
    spex_dynamic_to_CSC_mpz.o \
    spex_CSC_mpz_to_dynamic.o \
    spex_delete_mpz_array.o \
    SPEX_determine_symmetry.o SPEX_transpose.o \

$(CS_Util): ../SPEX_Util/Include/SPEX_Util.h ../SPEX_Util/Source/spex_util_internal.h

.PRECIOUS: SPEX_%.c spex_%.c

SPEX_%.c:
	- ln -s ../SPEX_Util/Source/$@

spex_%.c:
	- ln -s ../SPEX_Util/Source/$@

CS_Left_LU = \
    SPEX_LU_analyze.o \
    SPEX_Left_LU_backslash.o \
    spex_left_lu_back_sub.o \
    spex_left_lu_dfs.o \
    SPEX_Left_LU_factorize.o \
    spex_left_lu_forward_sub.o \
    spex_left_lu_get_largest_pivot.o \
    spex_left_lu_get_nonzero_pivot.o \
    spex_left_lu_get_pivot.o \
    spex_left_lu_get_smallest_pivot.o \
    spex_left_lu_permute_b.o \
    spex_left_lu_permute_x.o \
    spex_left_lu_reach.o \
    spex_left_lu_ref_triangular_solve.o \
    SPEX_Left_LU_solve.o

$(CS_Left_LU): ../SPEX_Left_LU/Include/SPEX_Left_LU.h ../SPEX_Left_LU/Source/spex_left_lu_internal.h ../SPEX_Left_LU/Demo/demos.h

.PRECIOUS: SPEX_Left_LU_%.c spex_left_lu_%.c spexlu_demo.c demos.c

SPEX_Left_LU_%.c:
	- ln -s ../SPEX_Left_LU/Source/$@

spex_left_lu_%.c:
	- ln -s ../SPEX_Left_LU/Source/$@

spexlu_demo.c:
	- ln -s ../SPEX_Left_LU/Demo/$@

demos.c:
	- ln -s ../SPEX_Left_LU/Demo/$@

test4llu: $(CS_Tcov) $(CS_Util) $(CS_Left_LU) tcov_for_left_LU.c demos.c
	$(CC) $(LDFLAGS) tcov_for_left_LU.c demos.c $(CFLAGS) -o tcov_for_left_LU $(CS_Tcov) $(CS_Util) $(CS_Left_LU) $(LDLIBS)

SPEXLLU: $(CS_Util) $(CS_Left_LU) spexlu_demo.c demos.c
	$(CC) $(LDFLAGS) spexlu_demo.c demos.c $(CFLAGS) -o spexlu_demo $(CS_Tcov) $(CS_Util) $(CS_Left_LU) $(LDLIBS)

CS_Update = \
    spex_update_get_scattered_v.o \
    spex_update_cppu.o \
    spex_update_dppu1.o \
    spex_update_dppu2.o \
    spex_update_insert_new_entry.o \
    spex_update_ipge.o \
    spex_update_finalize_and_insert_vk.o \
    spex_update_triangular_solve.o \
    spex_update_backward_sub.o \
    spex_update_forward_sub.o \
    spex_update_verify.o \
    SPEX_Update_LU_ColRep.o \
    SPEX_Update_Solve.o \
    SPEX_Update_permute_row.o \
    SPEX_Update_matrix_canonicalize.o \
    SPEX_Update_Chol_Rank1.o \


$(CS_Update): ../SPEX_Update/Include/SPEX_Update.h ../SPEX_Update//Source/spex_update_internal.h

.PRECIOUS: SPEX_Update_%.c spex_update_%.c LU_update_demo.c

SPEX_Update_%.c:
	- ln -s ../SPEX_Update/Source/$@

spex_update_%.c:
	- ln -s ../SPEX_Update/Source/$@

LU_update_demo.c:
	- ln -s ../SPEX_Update/Demo/$@

simple_rand.c:
	- ln -s ../SPEX_Update/Demo/$@

test4luu: $(CS_Tcov) $(CS_Util) $(CS_Update) tcov_for_update.c simple_rand.c
	$(CC) $(LDFLAGS) tcov_for_update.c simple_rand.c $(CFLAGS) -o tcov_for_update $(CS_Tcov) $(CS_Util) $(CS_Update) $(LDLIBS)

SPEXLUU: $(CS_Util) $(CS_Update) LU_update_demo.c
	$(CC) $(LDFLAGS) LU_update_demo.c $(CFLAGS) -o LU_update_demo $(CS_Util) $(CS_Update) $(LDLIBS)

# run all statement coverage tests
runtests: runtests4llu runtests4luu

# only run test for SPEX_Left_LU
runtests4llu: all
	- ./spexlu_demo p 2 q 0 f ../SPEX_Left_LU/ExampleMats/10teams_mat.txt ../SPEX_Left_LU/ExampleMats/10teams_v.txt
	- ./spexlu_demo p 3 q 1 o 1 f ../SPEX_Left_LU/ExampleMats/10teams_mat.txt ../SPEX_Left_LU/ExampleMats/10teams_v.txt
	- ./spexlu_demo p 4 q 2 o 1 f ../SPEX_Left_LU/ExampleMats/test_mat.txt ../SPEX_Left_LU/ExampleMats/test_rhs.txt
	- ./spexlu_demo p 5 f ../SPEX_Left_LU/ExampleMats/10teams_mat.txt ../SPEX_Left_LU/ExampleMats/10teams_v.txt
	- ./tcov_for_left_LU > test_for_left_LU.out
	- ./tcov_for_left_LU 0 1 1

# only run test for SPEX_Update
runtests4luu: all
	- ./tcov_for_update 0 6 > test_for_update1.out
	- ./tcov_for_update 0 0 25 > test_for_update2.out
	- ./tcov_for_update 1 6 > test_for_update3.out
	- #./LU_update_demo

# To run with valgrind:
V = valgrind #--leak-check=full

# run all statement coverage tests but with valgrind
vtests: vtests4llu vtests4luu

# run test for SPEX_Left_LU with valgrind
vtests4llu: all
	- $(V) ./spexlu_demo p 2 q 0 f ../SPEX_Left_LU/ExampleMats/10teams_mat.txt ../SPEX_Left_LU/ExampleMats/10teams_v.txt
	- $(V) ./spexlu_demo p 3 q 1 o 1 f ../SPEX_Left_LU/ExampleMats/10teams_mat.txt ../SPEX_Left_LU/ExampleMats/10teams_v.txt
	- $(V) ./spexlu_demo p 4 q 2 o 1 f ../SPEX_Left_LU/ExampleMats/test_mat.txt ../SPEX_Left_LU/ExampleMats/test_rhs.txt
	- $(V) ./spexlu_demo p 5 f ../SPEX_Left_LU/ExampleMats/10teams_mat.txt ../SPEX_Left_LU/ExampleMats/10teams_v.txt
	- $(V) ./tcov_for_left_LU > test_for_left_LU.out
	- $(V) ./tcov_for_left_LU 0 1 1

# run test for SPEX_Update with valgrind
vtests4luu: all
	- $(V) ./tcov_for_update 0 6 > test_for_update1.out
	- $(V) ./tcov_for_update 0 0 25 > test_for_update2.out
	- $(V) ./tcov_for_update 1 6 > test_for_update3.out

# remove all files not in the original distribution, including symbolic links
clean:
	- $(RM) *.o *.bbg *.da *.gcov *.gcda *gcno
	- $(RM) tcov_for_left_LU tcov_for_update spexlu_demo LU_update_demo *.out *.a cov.sort out Mats4Tcov/SPEX_Update/mat.txt
	- $(RM) -r SPEX_*.c spex_*.c *.dSYM spexlu_demo.c demos.c LU_update_demo.c simple_rand.c

purge: distclean

distclean: clean

