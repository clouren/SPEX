#-------------------------------------------------------------------------------
# SPEX_CHOLMOD/Tcov/Makefile: compile and run SPEX_CHOLMOD test coverage
#-------------------------------------------------------------------------------

# SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
# Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
# SPEX_CHOLMOD/License for the license.

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
        -I../../../SuiteSparse_config \
        -DSPEX_GMP_LIST_INIT=2 -DSPEX_LU_UPDATE_TCOV -DSPEX_DEBUG

LDLIBS += -lm -lgmp -lmpfr -lsuitesparseconfig

# run all statement coverage tests, and then check for 100% coverage
run: runtests
	./covall

all: test demo

CS2 = SPEX_initialize.o \
    SPEX_initialize_expert.o \
    SPEX_finalize.o \
    SPEX_realloc.o \
    SPEX_malloc.o \
    SPEX_calloc.o \
    SPEX_free.o \
    SPEX_gmp.o \
    SPEX_cumsum.o \
    SPEX_create_default_options.o \
    spex_create_mpq_array.o \
    spex_create_mpz_array.o  \

$(CS2): ../../SPEX_Util/Include/SPEX_Util.h ../../SPEX_Util/Source/spex_util_internal.h

.PRECIOUS: SPEX_%.c spex_%.c

SPEX_%.c:
	- ln -s ../../SPEX_Util/Source/$@

spex_%.c:
	- ln -s ../../SPEX_Util/Source/$@

CS1 = tcov_malloc_test.o \
    spex_get_scattered_v.o \
    spex_cppu.o \
    spex_dppu1.o \
    spex_dppu2.o \
    spex_insert_new_entry.o \
    spex_ipge.o \
    SPEX_LUU.o \
    spex_delete_mpq_array.o \
    spex_delete_mpz_array.o \
    SPEX_mat_alloc.o \
    SPEX_mat_free.o \
    spex_finalize_and_insert_vk.o \
    spex_triangular_solve.o \
    spex_find_next_nz.o \
    SPEX_vector_alloc.o \
    SPEX_vector_free.o \
    SPEX_vector_realloc.o \
    SPEX_solve.o \
    spex_backward_sub.o \
    spex_forward_sub.o \
    spex_get_nnz_pattern.o \
    spex_verify.o \
    SPEX_create_mpq_array.o \
    SPEX_create_mpz_array.o  \
    SPEX_delete_mpq_array.o \
    SPEX_delete_mpz_array.o  \
#spex_delete_scattered_vector.o  spex_create_scattered_vector.o

$(CS1): ../Include/SPEX_LU_Update.h ../Source/spex_lu_update_internal.h tcov_malloc_test.h #../../SPEX_Util/Include/SPEX_Util.h ../../SPEX_Util/Source/spex_util_internal.h

.PRECIOUS: SPEX_%.c spex_%.c

SPEX_%.c:
	- ln -s ../Source/$@

spex_%.c:
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

