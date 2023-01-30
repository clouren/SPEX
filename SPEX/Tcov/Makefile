#-------------------------------------------------------------------------------
# SPEX/Tcov/Makefile: compile and run SPEX test coverage
#-------------------------------------------------------------------------------

# SPEX: (c) 2019-2023, Chris Lourenco, Jinhao Chen,
# Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
# All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#-------------------------------------------------------------------------------

# FIXME: use the right demos

default: cov

SUITESPARSE ?= $(realpath $(CURDIR)/../..)

################################################################################

    INSTALL_LIB ?= $(SUITESPARSE)/lib

    # use 8 jobs by default
    JOBS ?= 8

    LDFLAGS = --coverage
    LDFLAGS += -L$(INSTALL_LIB)
    LDFLAGS += -Wl,-rpath=$(INSTALL_LIB)

    # remove object files, but keep compiled libraries via 'make clean'
    CLEAN = *.o *.obj *.ln *.bb *.bbg *.da *.tcov *.gcov gmon.out *.bak *.d \
        *.gcda *.gcno *.aux *.bbl *.blg *.log *.toc *.dvi *.lof *.lot

    # also remove compiled libraries, via 'make distclean'
    PURGE = *.so* *.a *.dll *.dylib *.dSYM

################################################################################

# select the thread-safety mechanism to test: none, openmp, or pthreads.
# THREADS =
# THREADS = -fopenmp
  THREADS = -DSPEX_USE_PTHREADS -pthread

# Linux test coverage (gcc is required for test coverage)
CC = gcc
CFLAGS = -g -fprofile-arcs -ftest-coverage \
        -Wall -W -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
        -Wredundant-decls -Wnested-externs -Wdisabled-optimization -std=c99 \
        -Wno-unused-parameter \
	-I../Include/ \
	-I../Demo/ \
	-I../SPEX_Utilities/Source/ \
	-I../SPEX_Cholesky/Source/ \
	-I../SPEX_LU/Source  \
	-I../SPEX_Update/Source \
        -I../../SuiteSparse_config -I../../COLAMD/Include \
	-I../../AMD/Include \
        -DSPEX_GMP_LIST_INIT=2 -DSPEX_GMP_TEST_COVERAGE $(THREADS)

LDLIBS += -lgmp -lmpfr -lcolamd -lamd -lsuitesparseconfig -lm

# run all statement coverage tests, and then check for 100% coverage
cov: runtests
	./covall

all: tcov_for_lu spexlu_demo tcov_for_update tcov_for_cholesky

#-------------------------------------------------------------------------------
# compile .c file in this folder

OBJ_Tcov = tcov_utilities.o

$(OBJ_Tcov): tcov_utilities.h

#-------------------------------------------------------------------------------
# link and compile each .c file in SPEX_Utilities/Source

OBJ_Util = \
    SPEX_version.o \
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
    spex_colamd.o \
    spex_create_mpq.o

$(OBJ_Util): ../Include/SPEX.h ../SPEX_Utilities/Source/spex_util_internal.h

.PRECIOUS: SPEX_%.c spex_%.c

SPEX_%.c:
	- ln -s ../SPEX_Utilities/Source/$@

spex_%.c:
	- ln -s ../SPEX_Utilities/Source/$@

#-------------------------------------------------------------------------------
# link and compile each .c file in SPEX_LU/Source

OBJ_LU = \
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

.PRECIOUS: SPEX_lu_%.c SPEX_lu_%.c spex_left_lu_%.c spexlu_demo.c lu_demos.c

SPEX_lu_%.c:
	- ln -s ../SPEX_LU/Source/$@

SPEX_lu_%.c:
	- ln -s ../SPEX_LU/Source/$@

spex_left_lu_%.c:
	- ln -s ../SPEX_LU/Source/$@

spexlu_demo.c:
	- ln -s ../SPEX_LU/Demo/$@

lu_demos.c:
	- ln -s ../SPEX_LU/Demo/$@

#-------------------------------------------------------------------------------
# link and compile each .c file in SPEX_Update/Source

OBJ_Update = \
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
    spex_update_verify.o


$(OBJ_Update): ../Include/SPEX.h ../SPEX_Update/Source/spex_update_internal.h

.PRECIOUS: SPEX_update_%.c spex_update_%.c

SPEX_update_%.c:
	- ln -s ../SPEX_Update/Source/$@

spex_update_%.c:
	- ln -s ../SPEX_Update/Source/$@

#-------------------------------------------------------------------------------
# link and compile each .c file in SPEX_Cholesky/Source

OBJ_Cholesky = \
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
    spex_cholesky_left_triangular_solve.o \
    spex_cholesky_permute_A.o \
    spex_cholesky_post.o \
    spex_cholesky_pre_left_factor.o \
    spex_cholesky_preorder.o \
    SPEX_cholesky_solve.o \
    spex_cholesky_symbolic_analysis.o \
    spex_cholesky_tdfs.o \
    spex_cholesky_up_factor.o \
    spex_cholesky_up_triangular_solve.o

$(OBJ_Cholesky): ../Include/SPEX.h

.PRECIOUS: SPEX_cholesky_%.c spex_cholesky_%.c

SPEX_cholesky_%.c:
	- ln -s ../SPEX_Cholesky/Source/$@

spex_cholesky_%.c:
	- ln -s ../SPEX_Cholesky/Source/$@

chol_demos.c:
	- ln -s ../SPEX_Cholesky/Demo/$@

#-------------------------------------------------------------------------------

tcov_for_lu: $(OBJ_Tcov) $(OBJ_Util) $(OBJ_LU) $(OBJ_Update) tcov_for_lu.c
	$(CC) $(LDFLAGS) tcov_for_lu.c $(CFLAGS) -o tcov_for_lu $(OBJ_Tcov) $(OBJ_Util) $(OBJ_LU) $(OBJ_Update) $(LDLIBS)

spexlu_demo: $(OBJ_Util) $(OBJ_LU) $(OBJ_Update) spexlu_demo.c lu_demos.c ../SPEX_LU/Demo/lu_demos.h
	$(CC) $(LDFLAGS) spexlu_demo.c lu_demos.c $(CFLAGS) -o spexlu_demo $(OBJ_Tcov) $(OBJ_Util) $(OBJ_LU) $(OBJ_Update) $(LDLIBS)

tcov_for_update: $(OBJ_Tcov) $(OBJ_Util) $(OBJ_LU) $(OBJ_Update) tcov_for_update.c simple_rand.c simple_rand.h
	$(CC) $(LDFLAGS) tcov_for_update.c simple_rand.c $(CFLAGS) -o tcov_for_update $(OBJ_Tcov) $(OBJ_Util) $(OBJ_LU) $(OBJ_Update) $(LDLIBS)

tcov_for_cholesky: $(OBJ_Tcov) $(OBJ_Util) $(OBJ_Cholesky) tcov_for_cholesky.c simple_rand.c simple_rand.h chol_demos.c
	$(CC) $(LDFLAGS) $(CFLAGS) tcov_for_cholesky.c simple_rand.c chol_demos.c -o tcov_for_cholesky $(OBJ_Tcov) $(OBJ_Util) $(OBJ_Cholesky) $(LDLIBS)

#-------------------------------------------------------------------------------

# run all statement coverage tests
runtests: runtests4llu runtests4luu runtests4chol

# only run test for SPEX_LU
runtests4llu: all
	- ./spex_demo_lu_extended p 2 q 1 f ../ExampleMats/10teams.mat.txt ../ExampleMats/10teams.rhs.txt > lu1.out
	- ./spex_demo_lu_extended p 3 q 2 o 1 f ../ExampleMats/10teams.mat.txt ../ExampleMats/10teams.rhs.txt > lu2.out
	- ./spex_demo_lu_extended p 4 q 3 o 1 f ../ExampleMats/test.mat.txt ../ExampleMats/test.rhs.txt > lu3.out
	- ./spex_demo_lu_extended p 5 f ../ExampleMats/10teams.mat.txt ../ExampleMats/10teams.rhs.txt > lu4.out
	- ./tcov_for_lu > test_for_lu.out
	- ./tcov_for_lu 0 1 1 > lu5.out

# only run test for SPEX_Update
runtests4luu: all
	- ./tcov_for_update 0 6 > test_for_update1.out
	- ./tcov_for_update 0 0 26 > test_for_update2.out
	- ./tcov_for_update 1 6 > test_for_update3.out

# only run test for SPEX_Cholesky
runtests4chol: all
	- ./tcov_for_cholesky ../ExampleMats/mesh1e1.mat.txt ../ExampleMats/mesh1e1.rhs.txt > test_for_cholesky.out

#-------------------------------------------------------------------------------

# To run with valgrind:
# V = valgrind
V = valgrind --leak-check=full --show-leak-kinds=all

# run all statement coverage tests but with valgrind
vtests: vtests4llu vtests4luu vtests4chol

# run test for SPEX_LU with valgrind
vtests4llu: all
	- $(V) ./spexlu_demo p 2 q 0 f ../ExampleMats/10teams_mat.txt ../ExampleMats/10teams_v.txt > lu1.out
	- $(V) ./spexlu_demo p 3 q 1 o 1 f ../ExampleMats/10teams_mat.txt ../ExampleMats/10teams_v.txt > lu2.out
	- $(V) ./spexlu_demo p 4 q 2 o 1 f ../ExampleMats/test_mat.txt ../ExampleMats/test_rhs.txt > lu3.out
	- $(V) ./spexlu_demo p 5 f ../ExampleMats/10teams_mat.txt ../ExampleMats/10teams_v.txt > lu4.out
	- $(V) ./tcov_for_lu > test_for_lu.out
	- $(V) ./tcov_for_lu 0 1 1 > lu5.out

# run test for SPEX_Update with valgrind
vtests4luu: all
	- $(V) ./tcov_for_update 0 6 > test_for_update1.out
	- $(V) ./tcov_for_update 0 0 25 > test_for_update2.out
	- $(V) ./tcov_for_update 1 6 > test_for_update3.out

# run test for SPEX_Cholesky with valgrind
vtests4chol: all
	- $(V) ./tcov_for_cholesky ../ExampleMats/872.mat.txt > test_for_cholesky.out

#-------------------------------------------------------------------------------

# remove all files not in the original distribution, including symbolic links
clean:
	- $(RM) *.o *.bbg *.da *.gcov *.gcda *gcno
	- $(RM) tcov_for_lu tcov_for_update spexlu_demo *.out *.a out Mats4Tcov/SPEX_Update/mat.txt
	- $(RM) -r SPEX_*.c spex_*.c *.dSYM spexlu_demo.c lu_demos.c
	- $(RM) -r tcov_for_cholesky chol_demos.c

purge: distclean

distclean: clean
