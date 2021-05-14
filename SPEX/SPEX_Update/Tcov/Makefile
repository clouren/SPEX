#-------------------------------------------------------------------------------
# SPEX_CHOLMOD/Tcov/Makefile: compile and run SPEX_CHOLMOD test coverage
#-------------------------------------------------------------------------------

# SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
# Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
# SPEX_CHOLMOD/License for the license.

#-------------------------------------------------------------------------------

default: run

SUITESPARSE ?= $(realpath $(CURDIR)/../../..)
include ../../../SuiteSparse_config/SuiteSparse_config.mk

# Linux test coverage
CC = gcc
CFLAGS = -g -fprofile-arcs -ftest-coverage \
        -Wall -W -Wshadow -Wmissing-prototypes -Wstrict-prototypes \
        -Wredundant-decls -Wnested-externs -Wdisabled-optimization -std=c99 \
        -Wno-unused-parameter -I../Include -I../Source -I../Demo \
        -I../../SPEX_Util/Include/ -I../../SPEX_Util/Source/ \
        -I../../../SuiteSparse_config \
        -DSPEX_GMP_LIST_INIT=2 -DSPEX_UPDATE_TCOV #-DSPEX_DEBUG

LDLIBS += -lm -lgmp -lmpfr -lsuitesparseconfig

# run all statement coverage tests, and then check for 100% coverage
run: runtests
	./covall

all: test demo

CS2 = \
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
    SPEX_malloc.o \
    SPEX_matrix_allocate.o \
    SPEX_matrix_check.o \
    SPEX_matrix_copy.o \
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
    spex_delete_mpz_array.o

$(CS2): ../../SPEX_Util/Include/SPEX_Util.h ../../SPEX_Util/Source/spex_util_internal.h

.PRECIOUS: SPEX_%.c spex_%.c

spex_%.c:
	- ln -s ../../SPEX_Util/Source/$@

SPEX_%.c:
	- ln -s ../../SPEX_Util/Source/$@

CS1 = tcov_malloc_test.o \
    spex_update_get_scattered_v.o \
    spex_update_cppu.o \
    spex_update_dppu1.o \
    spex_update_dppu2.o \
    spex_update_insert_new_entry.o \
    spex_update_ipge.o \
    SPEX_Update_LU_ColRep.o \
    spex_update_finalize_and_insert_vk.o \
    spex_update_triangular_solve.o \
    spex_update_find_next_nz.o \
    SPEX_Update_Solve.o \
    spex_update_backward_sub.o \
    spex_update_forward_sub.o \
    spex_update_verify.o \
#spex_delete_scattered_vector.o  spex_create_scattered_vector.o

$(CS1): ../Include/SPEX_Update.h ../Source/spex_update_internal.h tcov_malloc_test.h #../../SPEX_Util/Include/SPEX_Util.h ../../SPEX_Util/Source/spex_util_internal.h

.PRECIOUS: SPEX_%.c spex_%.c

SPEX_Update%.c:
	- ln -s ../Source/$@

spex_update%.c:
	- ln -s ../Source/$@

LU_update_demo.c:
	- ln -s ../Demo/$@

#demos.c:
#	- ln -s ../Demo/$@

test: $(CS1) $(CS2) tcov_test.c
	$(CC) $(LDFLAGS) tcov_test.c $(CFLAGS) -o tcov_test $(CS1) $(CS2) $(LDLIBS)

demo: $(CS1) $(CS2) LU_update_demo.c
	$(CC) $(LDFLAGS) LU_update_demo.c $(CFLAGS) -o LU_update_demo $(CS1) $(CS2) $(LDLIBS)

rundemo: demo
	- ./LU_update_demo > demo.out

# run all statement coverage tests
runtests: all
	- ./tcov_test 4 30 > test.out
	- ./tcov_test 0 11 > test.out
	- #./LU_update_demo

# To run with valgrind:
V = valgrind #--leak-check=full

# run all statement coverage tests but with valgrind
vtests: all
	- $(V) ./tcov_test
	- $(V) ./LU_update_demo

# remove all files not in the original distribution, including symbolic links
clean:
	- $(RM) *.o *.bbg *.da *.gcov *.gcda *gcno
	- $(RM) tcov_test LU_update_demo *.out *.a cov.sort out
	- $(RM) -r SPEX_*.c spex_*.c *.dSYM LU_update_demo.c

purge: distclean

distclean: clean

