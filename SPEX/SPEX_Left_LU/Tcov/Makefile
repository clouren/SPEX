#-------------------------------------------------------------------------------
# SPEX/SPEX/SPEX_Left_LU/Tcov/Makefile: compile and run SLIP_LU test coverage
#-------------------------------------------------------------------------------

# SPEX: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
# Timothy A. Davis, Texas A&M University.  All Rights Reserved.  See
# SPEX/License for the license.

#-------------------------------------------------------------------------------

default: run

include ../../../SuiteSparse_config/SuiteSparse_config.mk

# Linux test coverage
CC = gcc
CFLAGS = -g -fprofile-arcs -ftest-coverage \
        -Wall -W -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
        -Wredundant-decls -Wnested-externs -Wdisabled-optimization -std=c99 \
        -Wno-unused-parameter -I../Include -I../Source -I../Demo \
	-I../../SPEX_Util/Include/ -I../../SPEX_Util/Source/ \
        -I../../../SuiteSparse_config -I../../COLAMD/Include \
	-I../../AMD/Include \
        -DSPEX_GMP_LIST_INIT=2

LDLIBS += -lm -lgmp -lmpfr -lcolamd -lamd -lsuitesparseconfig

# run all statement coverage tests, and then check for 100% coverage
run: runtests
	./covall

all: test SPEXLLU

CS = \
    SPEX_calloc.o \
    SPEX_matrix_div.o \
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
    SPEX_LU_analysis_free.o \
    SPEX_malloc.o \
    SPEX_matrix_allocate.o \
    SPEX_matrix_check.o \
    SPEX_matrix_copy.o \
    SPEX_matrix_free.o \
    SPEX_matrix_mul.o \
    SPEX_matrix_nnz.o \
    SPEX_realloc.o \
    SPEX_sparse_collapse.o \
    SPEX_sparse_realloc.o \
    tcov_malloc_test.o \
    SPEX_Left_LU_analyze.o \
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
    SPEX_Left_LU_solve.o #SPEX_determine_symmetry.o SPEX_transpose.o \

$(CS): ../Include/SPEX_Left_LU.h ../../SPEX_Util/Include/SPEX_Util.h ../Source/spex_left_lu_internal.h ../../SPEX_Util/Source/spex_util_internal.h tcov_malloc_test.h ../Demo/demos.h

.PRECIOUS: SPEX_Left_LU_%.c spex_left_lu_%.c SPEX_Left_LU.c demos.c SPEX_%.c spex_%.c

SPEX_Left_LU_%.c:
	- ln -s ../Source/$@

spex_left_lu_%.c:
	- ln -s ../Source/$@

SPEX_Left_LU.c:
	- ln -s ../Demo/$@

demos.c:
	- ln -s ../Demo/$@

SPEX_%.c:
	- ln -s ../../SPEX_Util/Source/$@

spex_%.c:
	- ln -s ../../SPEX_Util/Source/$@

test: $(CS) tcov_test.c demos.c
	$(CC) $(LDFLAGS) tcov_test.c demos.c $(CFLAGS) -o tcov_test $(CS) $(LDLIBS)

SPEXLLU: $(CS) SPEX_Left_LU.c demos.c
	$(CC) $(LDFLAGS) SPEX_Left_LU.c demos.c $(CFLAGS) -o SPEX_Left_LU $(CS) $(LDLIBS)

# run all statement coverage tests
runtests: all
	- ./SPEX_Left_LU p 2 q 0
	- ./SPEX_Left_LU p 3 q 1 o 1
	- ./SPEX_Left_LU p 4 q 2 o 1 f ../ExampleMats/test_mat.txt ../ExampleMats/test_rhs.txt
	- ./SPEX_Left_LU p 5
	- ./tcov_test
	- ./tcov_test 0 1 1

# To run with valgrind:
V = valgrind #--leak-check=full

# run all statement coverage tests but with valgrind
vtests: all
	- $(V) ./SPEX_Left_LU p 2 q 0
	- $(V) ./SPEX_Left_LU p 3 q 1 o 1
	- $(V) ./SPEX_Left_LU p 4 q 2 o 1 f ../ExampleMats/test_mat.txt ../ExampleMats/test_rhs.txt
	- $(V) ./SPEX_Left_LU p 5
	- $(V) ./tcov_test
	- $(V) ./tcov_test 0 1 1

# remove all files not in the original distribution, including symbolic links
clean:
	- $(RM) *.o *.bbg *.da *.gcov *.gcda *gcno
	- $(RM) tcov_test SPEX_Left_LU *.out *.a cov.sort out
	- $(RM) -r SPEX_*.c spex_*.c *.dSYM SPEX_Left_LU.c demos.c

purge: distclean

distclean: clean

