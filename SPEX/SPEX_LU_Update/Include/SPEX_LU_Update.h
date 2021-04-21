//------------------------------------------------------------------------------
// SPEX_LU_Update/Include/SPEX_LU_Update.h: user #include file for SPEX_LU_Update.
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

#ifndef SPEX_LU_UPDATE_H
#define SPEX_LU_UPDATE_H

// This software package exactly updates the Cholesky factorization of a
// sparse matrix. This code accompanies the paper


//    The theory associated with this software can be found in the paper


//    If you use this code, you must first download and install the GMP and
//    MPFR libraries. GMP and MPFR can be found at:
//              https://gmplib.org/
//              http://www.mpfr.org/

//    If you use SPEX CHOLMOD for a publication, we request that you
//    please cite the above two papers.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Jinhao Chen, Timothy Davis, and Erick Moreno-Centeno
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Contact Information----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Please contact Jinhao Chen (cjh10644@hotmail.com)
//    or Tim Davis (timdavis@aldenmath.com, DrTimothyAldenDavis@gmail.com,
//                  davis@tamu.edu)

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Copyright--------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    SPEX CHOLMOD is free software; you can redistribute it and/or modify
//     it under the terms of either:
//
//        * the GNU Lesser General Public License as published by the
//          Free Software Foundation; either version 3 of the License,
//          or (at your option) any later version.
//
//     or
//
//        * the GNU General Public License as published by the Free Software
//          Foundation; either version 2 of the License, or (at your option) any
//          later version.
//
//    or both in parallel, as here.
//
//    See license.txt for license info.
//
// This software is copyright by Jinhao Chen, Timothy A. Davis and Erick
// Moreno-Centeno. All Rights Reserved.
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------DISCLAIMER-----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// SPEX CHOLMOD is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//--------------------------Summary---------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    This software package update a REF Cholesky factorization A = LDL^T
//    exactly when A is changed with a row and column. The input matrices are
//    stored as either integers, double precision numbers, multiple precision
//    floating points (through the mpfr library) or as rational numbers (as a
//    collection of numerators and denominators using the GMP mpq_t data
//    structure). Appropriate routines within the code transform the input into
//    an integral matrix in compressed column form.


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------Include files required by SPEX CHOLMOD---------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "SuiteSparse_config.h"
#include "SPEX_Util.h"

//------------------------------------------------------------------------------
// Version
//------------------------------------------------------------------------------

// Current version of the code
#define SPEX_LU_UPDATE_VERSION "1.0.0"
#define SPEX_LU_UPDATE_VERSION_MAJOR 1
#define SPEX_LU_UPDATE_VERSION_MINOR 0
#define SPEX_LU_UPDATE_VERSION_SUB   0

//------------------------------------------------------------------------------
// Error codes
//------------------------------------------------------------------------------

// Most SPEX_LU_UPDATE functions return a code that indicates if it was successful
// or not. Otherwise the code returns a pointer to the object that was created
// or it returns void (in the case that an object was deleted)
/*
typedef enum
{
    SPEX_OK = 0,              // all is well
    SPEX_OUT_OF_MEMORY = -1,  // out of memory
    SPEX_SINGULAR = -2,       // the input matrix A is singular
    SPEX_INCORRECT_INPUT = -3,// one or more input arguments are incorrect
    SPEX_PANIC = -4           // SPEX_CHOLMOD used without proper initialization
}
SPEX_info ;*/

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Data Structures--------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
/*
// This struct serves as a global struct to define all user-selectable options.

typedef struct SPEX_options
{
    int print_level ;      // 0: print nothing, 1: just errors,
                           // 2: terse (basic stats),
                           // 3: all, with matrices and results
    int32_t prec ;         // Precision used to output file if MPFR is chosen
    mpfr_rnd_t round ;     // Type of MPFR rounding used
} SPEX_options ;

// Purpose: Create and return SPEX_options object with default parameters
// upon successful allocation, which are defined in SPEX_CHOLMOD_internal.h
// To free it, simply use SPEX_FREE (option).
SPEX_options* SPEX_create_default_options (void) ;*/

//------------------------------------------------------------------------------
// SPEX_matrix: a sparse CSC, sparse triplet, or dense matrix
//------------------------------------------------------------------------------
typedef struct
{
    int64_t nz;   // number of nonzeros. nz is meaningless for a dense vector
    int64_t nzmax;// size of array i and x, nz <= nzmax
    int64_t *i;   // array of size nzmax that contains the column/row indices
                  // of each nnz. For a dense vector, i == NULL
    mpz_t *x;     // array of size nzmax that contains the values of each nnz
} SPEX_vector;

SPEX_info SPEX_vector_alloc 
( 
    SPEX_vector **v_handle,         // vector to be allocated
    const int64_t nzmax,            // number of nnz entries in v
    const bool IsSparse              // indicate if the vector is sparse
);

SPEX_info SPEX_vector_realloc
(
    SPEX_vector* v, // the vector to be expanded
    const int64_t new_size// desired new size for v
);

void SPEX_vector_free
(
    SPEX_vector **v  // vector to be deleted
);

typedef struct
{
    int64_t m ;   // number of entries in each vector
    int64_t n ;   // number of vectors
    SPEX_vector **v;// array of size n, each entry of this array is a column/row
                  // vector.
    mpq_t scale ; // the scale factor applied to original entries to make
                  // them integers. That is, the original value can be obtained
                  // as v->x[*]/scale
} SPEX_mat;

SPEX_info SPEX_mat_alloc
(
    SPEX_mat **A_handle,      // matrix to be allocated
    const int64_t n,             // number of vectors
    const int64_t m,             // number of max entries in each vector
    const bool IsSparse          // indicate if the vector is sparse
);

void SPEX_mat_free
(
    SPEX_mat **A  // matrix to be deleted
);

SPEX_info SPEX_CSC_to_mat
(
    SPEX_mat **A_handle,          // converted SPEX_mat matrix
    const int64_t *perm,          // vector for permutation matrix, consider as
                                  // identity matrix if input as NULL
    const bool A_Is_ColWise,      // true if A->v[i] is the i-th col of A.
                                  // Otherwise, A->v[i] is the i-th row of A
    const SPEX_matrix *B,         // original matrix
    const SPEX_options *option
);

SPEX_info SPEX_mat_to_CSC
(
    SPEX_matrix **A_handle,       // converted CSC matrix
    const SPEX_mat *B,            // original matrix
    const int64_t *perm_inv,      // inverse of permutation, consider as
                                  // identity matrix if input as NULL
    const bool B_Is_ColWise,      // true if B->v[i] is the i-th col of B.
                                  // Otherwise, B->v[i] is the i-th row of B
    const SPEX_options *option
);

mpz_t* SPEX_create_mpz_array
(
    int64_t n            // size of the array
);

void SPEX_delete_mpz_array
(
    mpz_t **x,      // mpz array to be deleted
    int64_t n       // Size of x
);

mpq_t* SPEX_create_mpq_array
(
    int64_t n            // size of the array
);

void SPEX_delete_mpq_array
(
    mpq_t **x,      // mpq array to be deleted
    int64_t n       // Size of x
);
/*
//------------------------------------------------------------------------------
// SPEX_matrix: a sparse CSC, sparse triplet, or dense matrix
//------------------------------------------------------------------------------

// SPEX_CHOLMOD uses a single matrix data type, SPEX_matrix, which can be held in
// one of three kinds of formats:  sparse CSC (compressed sparse column),
// sparse triplet, and dense:

typedef enum
{
    SPEX_CSC = 0,           // matrix is in compressed sparse column format
    SPEX_CSR = 1,           // matrix is in compressed sparse row    format
    SPEX_TRIPLET = 2,       // matrix is in sparse triplet format
    SPEX_DENSE = 3          // matrix is in dense format
}
SPEX_kind ;

// Each of the three formats can have values of 5 different data types: mpz_t,
// mpq_t, mpfr_t, int64_t, and double:

typedef enum
{
    SPEX_MPZ = 0,           // matrix of mpz_t integers
    SPEX_MPQ = 1,           // matrix of mpq_t rational numbers
    SPEX_MPFR = 2,          // matrix of mpfr_t
    SPEX_INT64 = 3,         // matrix of int64_t integers
    SPEX_FP64 = 4           // matrix of doubles
}
SPEX_type ;

// This gives a total of 15 different matrix types.  Not all functions accept
// all 15 matrices types, however.

// Suppose A is an m-by-n matrix with nz <= nzmax entries.
// The p, i, j, and x components are defined as:

// (0) SPEX_CSC:  A sparse matrix in CSC (compressed sparse column) format.
//      A->p is an int64_t array of size n+1, A->i is an int64_t array of size
//      nzmax (with nz <= nzmax), and A->x.type is an array of size nzmax of
//      matrix entries ('type' is one of mpz, mpq, mpfr, int64, or fp64).  The
//      row indices of column j appear in A->i [A->p [j] ... A->p [j+1]-1], and
//      the values appear in the same locations in A->x.type.  The A->j array
//      is NULL.  A->nz is ignored; nz is A->p [A->n].

// (1) SPEX_TRIPLET:  A sparse matrix in triplet format.  A->i and A->j are
//      both int64_t arrays of size nzmax, and A->x.type is an array of values
//      of the same size.  The kth tuple has row index A->i [k], column index
//      A->j [k], and value A->x.type [k], with 0 <= k < A->nz.  The A->p array
//      is NULL.

// (2) SPEX_DENSE:  A dense matrix.  The integer arrays A->p, A->i, and A->j
//      are all NULL.  A->x.type is a pointer to an array of size m*n, stored
//      in column-oriented format.  The value of A(i,j) is A->x.type [p]
//      with p = i + j*A->m.  A->nz is ignored; nz is A->m * A->n.

// The SPEX_matrix may contain 'shallow' components, A->p, A->i, A->j, and
// A->x.  For example, if A->p_shallow is true, then a non-NULL A->p is a
// pointer to a read-only array, and the A->p array is not freed by
// SPEX_matrix_free.  If A->p is NULL (for a triplet or dense matrix), then
// A->p_shallow has no effect.

typedef struct
{
    int64_t m ;         // number of rows
    int64_t n ;         // number of columns
    int64_t nzmax ;     // size of A->i, A->j, and A->x
    int64_t nz ;        // # nonzeros in a triplet matrix .
                        // Ignored for CSC, CSR and dense matrices.
    SPEX_kind kind ;    // CSC, CSR, triplet, or dense
    SPEX_type type ;    // mpz, mpq, mpfr, int64, or fp64 (double)

    bool fixed_size;    // if true, the nnz pattern can be purely determined
                        // by p. Otherwise

    int64_t *pnz;       // if CSC or CSR: column/row nnz, an array size is n.
                        // pnz[i] <= pend[i]-p[i]+1
                        // if triplet or dense or not fixed_size: A->p is NULL.
    int64_t *pend;      // If CSC or CSR: an array of size n,
                        // p[i] and pend[i] define the first and last slots for
                        // entries in i-th column/row. If triplet or dense or
                        // not fixed_size: A->p is NULL.

    int64_t *p ;        // if CSC or CSR: column/row pointers, an array size
                        // is n+1. if triplet or dense: A->p is NULL.
                        // if fixed_size is false, p[n] is useless
    bool p_shallow ;    // if true, A->p is shallow.

    int64_t *i ;        // if CSC or triplet: row indices, of size nzmax.
                        // if CSR or dense: A->i is NULL.
    bool i_shallow ;    // if true, A->i is shallow.

    int64_t *j ;        // if CSR or triplet: column indices, of size nzmax.
                        // if CSC or dense: A->j is NULL.
    bool j_shallow ;    // if true, A->j is shallow.

    union               // A->x.type has size nzmax.
    {
        mpz_t *mpz ;            // A->x.mpz
        mpq_t *mpq ;            // A->x.mpq
        mpfr_t *mpfr ;          // A->x.mpfr
        int64_t *int64 ;        // A->x.int64
        double *fp64 ;          // A->x.fp64
    } x ;
    bool x_shallow ;    // if true, A->x.type is shallow.

    mpq_t scale ;       // scale factor for mpz matrices (never shallow)
                        // For all matrices who's type is not mpz,
                        // mpz_scale = 1. 

} SPEX_matrix ;

//------------------------------------------------------------------------------
// SPEX_matrix_allocate: allocate an m-by-n SPEX_matrix
//------------------------------------------------------------------------------

// if shallow is false: All components (p,i,j,x) are allocated and set to zero,
//                      and then shallow flags are all false.

// if shallow is true:  All components (p,i,j,x) are NULL, and their shallow
//                      flags are all true.  The user can then set A->p,
//                      A->i, A->j, and/or A->x accordingly, from their own
//                      arrays.

SPEX_info SPEX_matrix_allocate
(
    SPEX_matrix **A_handle, // matrix to allocate
    SPEX_kind kind,         // CSC, triplet, or dense
    SPEX_type type,         // mpz, mpq, mpfr, int64, or double
    int64_t m,              // # of rows
    int64_t n,              // # of columns
    int64_t nzmax,          // max # of entries
    bool shallow,           // if true, matrix is shallow.  A->p, A->i, A->j,
                            // A->x are all returned as NULL and must be set
                            // by the caller.  All A->*_shallow are returned
                            // as true.
    bool init,              // If true, and the data types are mpz, mpq, or
                            // mpfr, the entries are initialized (using the
                            // appropriate SPEX_mp*_init function). If false,
                            // the mpz, mpq, and mpfr arrays are allocated but
                            // not initialized.
    const SPEX_options *option
) ;

//------------------------------------------------------------------------------
// SPEX_matrix_free: free a SPEX_matrix
//------------------------------------------------------------------------------

SPEX_info SPEX_matrix_free
(
    SPEX_matrix **A_handle, // matrix to free
    const SPEX_options *option
) ;

//------------------------------------------------------------------------------
// SPEX_matrix_nnz: # of entries in a matrix
//------------------------------------------------------------------------------

int64_t SPEX_matrix_nnz     // return # of entries in A, or -1 on error
(
    const SPEX_matrix *A,         // matrix to query
    const SPEX_options *option
) ;

//------------------------------------------------------------------------------
// SPEX_matrix_copy: makes a copy of a matrix
//------------------------------------------------------------------------------

// SPEX_matrix_copy: make a copy of a SPEX_matrix, into another kind and type.

SPEX_info SPEX_matrix_copy
(
    SPEX_matrix **C_handle, // matrix to create (never shallow)
    // inputs, not modified:
    SPEX_kind C_kind,       // C->kind: CSC, triplet, or dense
    SPEX_type C_type,       // C->type: mpz_t, mpq_t, mpfr_t, int64_t, or double
    SPEX_matrix *A,         // matrix to make a copy of (may be shallow)
    const SPEX_options *option
) ;

//------------------------------------------------------------------------------
// SPEX_matrix macros
//------------------------------------------------------------------------------

// These macros simplify the access to entries in a SPEX_matrix.
// The type parameter is one of: mpq, mpz, mpfr, int64, or fp64.

// To access the kth entry in a SPEX_matrix using 1D linear addressing,
// in any matrix kind (CSC, triplet, or dense), in any type:
#define SPEX_1D(A,k,type) ((A)->x.type [k])

// To access the (i,j)th entry in a 2D SPEX_matrix, in any type:
#define SPEX_2D(A,i,j,type) SPEX_1D (A, (i)+(j)*((A)->m), type)
*/
#define SPEX_LUU_2D(A,i,j) (A[(i)+(j)*(3)-1])
//------------------------------------------------------------------------------
// Memory management
//------------------------------------------------------------------------------

// SPEX_CHOLMOD relies on the SuiteSparse memory management functions,
// SuiteSparse_malloc, SuiteSparse_calloc, SuiteSparse_realloc, and
// SuiteSparse_free.
/*
// Allocate and initialize memory space for SPEX_CHOLMOD.
void *SPEX_calloc
(
    size_t nitems,      // number of items to allocate
    size_t size         // size of each item
) ;

// Allocate memory space for SPEX_CHOLMOD.
void *SPEX_malloc
(
    size_t size        // size of memory space to allocate
) ;

// Free the memory allocated by SPEX_calloc, SPEX_malloc, or SPEX_realloc.
void SPEX_free
(
    void *p         // pointer to memory space to free
) ;

// Free a pointer and set it to NULL.
#define SPEX_FREE(p)                        \
{                                           \
    SPEX_free (p) ;                         \
    (p) = NULL ;                            \
}

// SPEX_realloc is a wrapper for realloc.  If p is non-NULL on input, it points
// to a previously allocated object of size old_size * size_of_item.  The
// object is reallocated to be of size new_size * size_of_item.  If p is NULL
// on input, then a new object of that size is allocated.  On success, a
// pointer to the new object is returned.  If the reallocation fails, p is not
// modified, and a flag is returned to indicate that the reallocation failed.
// If the size decreases or remains the same, then the method always succeeds
// (ok is returned as true).

// Typical usage:  the following code fragment allocates an array of 10 int's,
// and then increases the size of the array to 20 int's.  If the SPEX_malloc
// succeeds but the SPEX_realloc fails, then the array remains unmodified,
// of size 10.
//
//      int *p ;
//      p = SPEX_malloc (10 * sizeof (int)) ;
//      if (p == NULL) { error here ... }
//      printf ("p points to an array of size 10 * sizeof (int)\n") ;
//      bool ok ;
//      p = SPEX_realloc (20, 10, sizeof (int), p, &ok) ;
//      if (ok) printf ("p has size 20 * sizeof (int)\n") ;
//      else printf ("realloc failed; p still has size 10 * sizeof (int)\n") ;
//      SPEX_free (p) ;

void *SPEX_realloc      // pointer to reallocated block, or original block
                        // if the realloc failed
(
    int64_t nitems_new,     // new number of items in the object
    int64_t nitems_old,     // old number of items in the object
    size_t size_of_item,    // sizeof each item
    void *p,                // old object to reallocate
    bool *ok                // true if success, false on failure
) ;

//------------------------------------------------------------------------------
// SPEX CHOLMOD memory environment routines
//------------------------------------------------------------------------------

// SPEX_initialize: initializes the working evironment for SPEX CHOLMOD library.
// It must be called prior to calling any other SPEX_* function.
SPEX_info SPEX_initialize (void) ;

// SPEX_initialize_expert is the same as SPEX_initialize, except that it allows
// for a redefinition of custom memory functions that are used for SPEX CHOLMOD
// and GMP.  The four inputs to this function are pointers to four functions
// with the same signatures as the ANSI C malloc, calloc, realloc, and free.
SPEX_info SPEX_initialize_expert
(
    void* (*MyMalloc) (size_t),             // user-defined malloc
    void* (*MyCalloc) (size_t, size_t),     // user-defined calloc
    void* (*MyRealloc) (void *, size_t),    // user-defined realloc
    void  (*MyFree) (void *)                // user-defined free
) ;

// SPEX_finalize: This function finalizes the working evironment for SPEX CHOLMOD
// library, and frees any internal workspace created by SPEX CHOLMOD.  It must be
// called as the last SPEX_* function called.
SPEX_info SPEX_finalize (void) ;*/

//------------------------------------------------------------------------------
// Primary factorization update routine
//------------------------------------------------------------------------------

SPEX_info SPEX_LUU
(
    SPEX_mat *A,         // the original matrix in compressed-column form
    SPEX_mat *L,         // stored in compressed-column form
    SPEX_mat *U,         // stored in comptessed-row form
    mpz_t *d,               // an array of size n that stores the unscaled pivot
    mpz_t *sd,              // an array of size n that stores the scaled pivot
    mpq_t *S,               // an array of size 3*n that stores pending scales
    int64_t *P,             // row permutation
    int64_t *P_inv,         // inverse of row permutation
    int64_t *Q,             // column permutation
    int64_t *Q_inv,         // inverse of column permutation
    SPEX_vector **vk,       // pointer to the inserted column, which will be
                            // swapped with A->v[k] in the output if succeed
    int64_t k,              // the column index that vk will be inserted
    const SPEX_options *option// command parameters
);

//------------------------------------------------------------------------------
// Function for solving LDUx =b
//------------------------------------------------------------------------------

SPEX_info SPEX_solve     // solves the linear system LD^(-1)U x = b
(
    // Output
    SPEX_mat **x_handle, // rational solution to the system
    // input:
    SPEX_mat *b,         // right hand side vector
    int64_t *h,             // history vector
    const SPEX_mat *A,   // Input matrix
    const SPEX_mat *L,   // lower triangular matrix
    const SPEX_mat *U,   // upper triangular matrix
    mpq_t *S,               // the pending scale factor matrix
    const mpz_t *sd,        // array of scaled pivots
    mpz_t *d,               // array of unscaled pivots
    const int64_t *Ldiag,   // L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const int64_t *Ucp,     // col pointers for col-wise nnz pattern of U
    const int64_t *Ucx,     // the value of k-th entry is found as 
                            // U->v[Uci[k]]->x[Ucx[k]]
    const int64_t *P,       // row permutation
    const int64_t *P_inv,   // inverse of row permutation
    const int64_t *Q,       // column permutation
    const int64_t *Q_inv,   // inverse of column permutation
    const bool keep_b,      // indicate if b will be reused
    const SPEX_options* option // Command options
);


SPEX_info SPEX_get_nnz_pattern    // find the nnz pattern of L and U
(
    // OUTPUT:
    int64_t **Ldiag,              // L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    int64_t **Lr_offdiag,         // Lr_offdiag[k] gives the column index of the
                                  // last off-diagonal nnz in k-th row of L.
                                  // -1 if no off diagonal entry
    int64_t **Uci,                // the row index for col-wise nnz pattern of U
    int64_t **Ucp,                // col pointers for col-wise pattern of U
    int64_t **Ucx,                // find the value of k-th entry as
                                  // U->v[Uci[k]]->x[Ucx[k]]
    // INPUT:
    const SPEX_mat *L,         // the target matrix L
    const SPEX_mat *U,         // the target matrix U
    const int64_t *P,             // row permutation
    const SPEX_options *option     // command option
);
//------------------------------------------------------------------------------
//---------------------------SPEX GMP/MPFR Functions----------------------------
//------------------------------------------------------------------------------

// The following functions are the SPEX CHOLMOD interface to the GMP/MPFR libary.
// Each corresponding GMP/MPFR function is given a wrapper to ensure that no
// memory leaks or crashes occur. All covered GMP functions can be found in
// SPEX_gmp.c

// The GMP library does not handle out-of-memory failures.  However, it does
// provide a mechanism for passing function pointers that replace GMP's use of
// malloc, realloc, and free.  This mechanism is used to provide a try/catch
// mechanism for memory allocation errors, using setjmp and longjmp.

// When a GMP function is called, this wrapper keeps track of a list of objects
// allocated by that function.  The list is started fresh each time a GMP
// function is called.  If any allocation fails, the NULL pointer is not
// returned to GMP.  Instead, all allocated blocks in the list are freed,
// and spex_gmp_allocate returns directly to wrapper.
/*
SPEX_info SPEX_mpfr_asprintf (char **str, const char *format, ... ) ;

SPEX_info SPEX_gmp_fscanf (FILE *fp, const char *format, ... ) ;

SPEX_info SPEX_mpz_init (mpz_t x) ;

SPEX_info SPEX_mpz_init2(mpz_t x, const size_t size) ;

SPEX_info SPEX_mpz_set (mpz_t x, const mpz_t y) ;

SPEX_info SPEX_mpz_set_ui (mpz_t x, const uint64_t y) ;

SPEX_info SPEX_mpz_set_si (mpz_t x, const int64_t y) ;

SPEX_info SPEX_mpz_swap (mpz_t x, mpz_t y);

SPEX_info SPEX_mpz_get_d (double *x, const mpz_t y) ;

SPEX_info SPEX_mpz_get_si (int64_t *x, const mpz_t y) ;

SPEX_info SPEX_mpz_set_q (mpz_t x, const mpq_t y) ;

SPEX_info SPEX_mpz_sub (mpz_t a, const mpz_t b, const mpz_t c) ;

SPEX_info SPEX_mpz_add (mpz_t a, const mpz_t b, const mpz_t c) ;

SPEX_info SPEX_mpz_addmul (mpz_t x, const mpz_t y, const mpz_t z) ;

SPEX_info SPEX_mpz_mul (mpz_t a, const mpz_t b, const mpz_t c) ;

SPEX_info SPEX_mpz_submul (mpz_t x, const mpz_t y, const mpz_t z) ;

SPEX_info SPEX_mpz_fdiv_q (mpz_t q, const mpz_t n, const mpz_t d) ;

SPEX_info SPEX_mpz_cdiv_q (mpz_t q, const mpz_t n, const mpz_t d) ;

SPEX_info SPEX_mpz_divexact (mpz_t x, const mpz_t y, const mpz_t z) ;

SPEX_info SPEX_mpz_gcd (mpz_t x, const mpz_t y, const mpz_t z) ;

SPEX_info SPEX_mpz_lcm (mpz_t lcm, const mpz_t x, const mpz_t y) ;

SPEX_info SPEX_mpz_neg (mpz_t x, const mpz_t y) ;

SPEX_info SPEX_mpz_abs (mpz_t x, const mpz_t y) ;

SPEX_info SPEX_mpz_cmp (int *r, const mpz_t x, const mpz_t y) ;

SPEX_info SPEX_mpz_cmpabs (int *r, const mpz_t x, const mpz_t y) ;

SPEX_info SPEX_mpz_cmp_ui (int *r, const mpz_t x, const uint64_t y) ;

SPEX_info SPEX_mpz_cmpabs_ui (int *r, const mpz_t x, const uint64_t y) ;

SPEX_info SPEX_mpz_sgn (int *sgn, const mpz_t x) ;

SPEX_info SPEX_mpz_sizeinbase (size_t *size, const mpz_t x, int64_t base) ;

SPEX_info SPEX_mpq_init (mpq_t x) ;

SPEX_info SPEX_mpq_set (mpq_t x, const mpq_t y) ;

SPEX_info SPEX_mpq_set_z (mpq_t x, const mpz_t y) ;

SPEX_info SPEX_mpq_canonicalize (mpq_t x);

SPEX_info SPEX_mpq_set_d (mpq_t x, const double y) ;

SPEX_info SPEX_mpq_swap (mpq_t x, mpq_t y);

SPEX_info SPEX_mpq_set_ui (mpq_t x, const uint64_t y, const uint64_t z) ;

SPEX_info SPEX_mpq_set_si (mpq_t x, const int64_t y, const uint64_t z) ;

SPEX_info SPEX_mpq_set_num (mpq_t x, const mpz_t y) ;

SPEX_info SPEX_mpq_set_den (mpq_t x, const mpz_t y) ;

SPEX_info SPEX_mpq_get_den (mpz_t x, const mpq_t y) ;

SPEX_info SPEX_mpq_get_d (double *x, const mpq_t y) ;

SPEX_info SPEX_mpq_neg (mpq_t x, const mpq_t y) ;

SPEX_info SPEX_mpq_abs (mpq_t x, const mpq_t y) ;

SPEX_info SPEX_mpq_add (mpq_t x, const mpq_t y, const mpq_t z) ;

SPEX_info SPEX_mpq_mul (mpq_t x, const mpq_t y, const mpq_t z) ;

SPEX_info SPEX_mpq_div (mpq_t x, const mpq_t y, const mpq_t z) ;

SPEX_info SPEX_mpq_cmp (int *r, const mpq_t x, const mpq_t y) ;

SPEX_info SPEX_mpq_cmp_ui (int *r, const mpq_t x,
                    const uint64_t num, const uint64_t den) ;

SPEX_info SPEX_mpq_cmp_z (int *r, const mpq_t x, const mpz_t y) ;

SPEX_info SPEX_mpq_sgn (int *sgn, const mpq_t x) ;

SPEX_info SPEX_mpq_equal (int *r, const mpq_t x, const mpq_t y) ;

SPEX_info SPEX_mpfr_init2(mpfr_t x, const uint64_t size) ;

SPEX_info SPEX_mpfr_set (mpfr_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_set_d (mpfr_t x, const double y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_set_si (mpfr_t x, int64_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_set_q (mpfr_t x, const mpq_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_set_z (mpfr_t x, const mpz_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_get_z (mpz_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_get_q (mpq_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_get_d (double *x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_get_si (int64_t *x, const mpfr_t y, const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_mul (mpfr_t x, const mpfr_t y, const mpfr_t z,
                    const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_mul_d (mpfr_t x, const mpfr_t y, const double z,
                    const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_div_d (mpfr_t x, const mpfr_t y, const double z,
                    const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_ui_pow_ui (mpfr_t x, const uint64_t y, const uint64_t z,
                    const mpfr_rnd_t rnd) ;

SPEX_info SPEX_mpfr_sgn (int *sgn, const mpfr_t x) ;

SPEX_info SPEX_mpfr_free_cache (void) ;

SPEX_info SPEX_mpfr_free_str (char *str) ;

SPEX_info SPEX_gmp_printf (const char *format, ... ) ;*/
#if 0
// These functions are currently unused, but kept here for future reference.
SPEX_info SPEX_gmp_asprintf (char **str, const char *format, ... ) ;
SPEX_info SPEX_mpfr_printf ( const char *format, ... ) ;
SPEX_info SPEX_gmp_fprintf (FILE *fp, const char *format, ... ) ;
SPEX_info SPEX_mpfr_fprintf (FILE *fp, const char *format, ... ) ;
SPEX_info SPEX_mpz_set_d (mpz_t x, const double y) ;
SPEX_info SPEX_mpfr_log2(mpfr_t x, const mpfr_t y, const mpfr_rnd_t rnd) ;
#endif

#endif

