//------------------------------------------------------------------------------
// SPEX_WAMF/SPEX_WAMF.h: user #include file for SPEX_WAMF
//------------------------------------------------------------------------------

// SPEX_WAMF: (c) 2021, Chris Lourenco, United States Naval Academy, Erick 
//            Moreno-Centeno Texas A&M
// All Rights Reserved.  See SPEX_WAMF/License for the license.

//------------------------------------------------------------------------------

#ifndef SPEX_WAMF_H
#define SPEX_WAMF_H


// This software performs a weighted approximate minimum fill ordering
// WARNING: This code is experimental and developmental, please do not use it.


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Christopher Lourenco and Erick Moreno-Centeno

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Contact Information----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Please contact Chris Lourenco (chrisjlourenco@gmail.com, lourenco@usna.edu)
//    or Tim Davis (timdavis@aldenmath.com, DrTimothyAldenDavis@gmail.com,
//                  davis@tamu.edu)


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Copyright--------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    SPEX_WAMF is free software; you can redistribute it and/or modify
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
// This software is copyright by Christopher Lourenco
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------DISCLAIMER-----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// SPEX WAMF is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//--------------------------Summary---------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// This software package computes a weighted AMF ordering of a matrix A. This ordering 
// is used to preorder sparse linear systems prior to exact factorization 
//
// This code accompanies the paper (submitted to XXX)

//    "A Weighted Approximate Minimum Fill Algorithm for Reducing Work Required in 
//     Integer-Preserving Factorization".

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------Include files required by SPEX WAMF----------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "SPEX_Util.h"

// SuiteSparse headers
#include "SuiteSparse_config.h"
#include "colamd.h"
#include "amd.h"

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Default Parameters-----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Current version of the code
#define SPEX_WAMF_VERSION "0.0.1"
#define SPEX_WAMF_VERSION_MAJOR 0
#define SPEX_WAMF_VERSION_MINOR 0
#define SPEX_WAMF_VERSION_SUB   1

#define ASSERT assert

// Default value of alpha
#define SPEX_WAMF_DEFAULT_ALPHA 0.8

// Note that WAMF inherits the same error and ordering codes as SLIP LU, thus they are omitted here


#ifndef FREE_WORKSPACE
#define FREE_WORKSPACE  
#endif

#define SPEX_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SPEX_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define SPEX_FLIP(i) (-(i)-2)
#define SPEX_UNFLIP(i) (((i) < 0) ? SLIP_FLIP(i) : (i))
#define SPEX_MARKED(Ap,j) (Ap [j] < 0)
#define SPEX_MARK(Ap,j) { Ap [j] = SLIP_FLIP (Ap [j]) ; }


// Field access macros for MPZ/MPQ/MPFR struct
// (similar definition in gmp-impl.h and mpfr-impl.h)

#define MPZ_SIZ(x)   ((x)->_mp_size)
#define MPZ_PTR(x)   ((x)->_mp_d)
#define MPZ_ALLOC(x) ((x)->_mp_alloc)
#define MPQ_NUM(x)   mpq_numref(x)
#define MPQ_DEN(x)   mpq_denref(x)
#define MPFR_MANT(x) ((x)->_mpfr_d)
#define MPFR_EXP(x)  ((x)->_mpfr_exp)
#define MPFR_PREC(x) ((x)->_mpfr_prec)
#define MPFR_SIGN(x) ((x)->_mpfr_sign)
#define MPFR_REAL_PTR(x) (&((x)->_mpfr_d[-1])) /*re-define but same result*/
/* Invalid exponent value (to track bugs...) */
#define MPFR_EXP_INVALID \
 ((mpfr_exp_t) 1 << (GMP_NUMB_BITS*sizeof(mpfr_exp_t)/sizeof(mp_limb_t)-2))

/* Macros to set the pointer in mpz_t/mpq_t/mpfr_t variable to NULL. It is best
 * practice to call these macros immediately after mpz_t/mpq_t/mpfr_t variable
 * is declared, and before the mp*_init function is called. It would help to
 * prevent error when SLIP_MP*_CLEAR is called before the variable is
 * successfully initialized.
 */

#define SPEX_MPZ_SET_NULL(x)                \
    MPZ_PTR(x) = NULL;                      \
    MPZ_SIZ(x) = 0;                         \
    MPZ_ALLOC(x) = 0;                       \

#define SPEX_MPQ_SET_NULL(x)                \
    MPZ_PTR(MPQ_NUM(x)) = NULL;             \
    MPZ_SIZ(MPQ_NUM(x)) = 0;                \
    MPZ_ALLOC(MPQ_NUM(x)) = 0;              \
    MPZ_PTR(MPQ_DEN(x)) = NULL;             \
    MPZ_SIZ(MPQ_DEN(x)) = 0;                \
    MPZ_ALLOC(MPQ_DEN(x)) = 0;              \

#define SPEX_MPFR_SET_NULL(x)               \
    MPFR_MANT(x) = NULL;                    \
    MPFR_PREC(x) = 0;                       \
    MPFR_SIGN(x) = 1;                       \
    MPFR_EXP(x) = MPFR_EXP_INVALID;         \

/* GMP does not give a mechanism to tell a user when an mpz, mpq, or mpfr
 * item has been cleared; thus, if mp*_clear is called on an object that
 * has already been cleared, gmp will crash. It is also not possible to
 * set a mp*_t = NULL. Thus, this mechanism modifies the internal GMP
 * size of entries to avoid crashing in the case that a mp*_t is cleared
 * multiple times.
 */

#define SPEX_MPZ_CLEAR(x)                   \
{                                           \
    if ((x) != NULL && MPZ_PTR(x) != NULL)  \
    {                                       \
        mpz_clear(x);                       \
        SPEX_MPZ_SET_NULL(x);               \
    }                                       \
}                                           \

#define SPEX_MPQ_CLEAR(x)                   \
{                                           \
    SPEX_MPZ_CLEAR(MPQ_NUM(x));             \
    SPEX_MPZ_CLEAR(MPQ_DEN(x));             \
}                                           \

#define SPEX_MPFR_CLEAR(x)                  \
{                                           \
    if ((x) != NULL && MPFR_MANT(x) != NULL)\
    {                                       \
        mpfr_clear(x);                      \
        SPEX_MPFR_SET_NULL(x);              \
    }                                       \
}                                           \

#define OK(method)                      \
{                                       \
    ok = method ;                       \
    if (ok != SPEX_OK)                  \
    {                                   \
        FREE_WORKSPACE ;                \
        return 0 ;                      \
    }                                   \
}                                       \

#define SPEX_CHECK(method)              \
{                                       \
    ok = method ;                       \
    if (ok != SPEX_OK)                  \
    {                                   \
        return 0 ;                      \
    }                                   \
}                                       \






/* Purpose: Given a matrix A which is stored as a sparse CSC matrix with GMP entries, 
 * convert it to a matrix of bit-lengths. That is the output matrix, B, has the same nonzero
 * pattern as A, but the entries are B(i,j) = bit-length( A(i,j)).
 * 
 * On success, this function returns the new matrix B. On failure, this function returns
 * a NULL pointer.
 */

SPEX_matrix* SPEX_WAMF_get_bitmat
(
    SPEX_matrix *A // matrix to be deleted
);


/* Purpose: WAMF stores a weight associated with each potential pivot row/column
 * (similar to the degree in minimum degree). This function initializes the weights
 * by creating an initial set of weights of the graph. 
 * 
 * Input arguments:
 * 
 *  A: matrix to be analyzed
 * 
 *  option: an int that determines how the weights are calculated. Currently these are:
 *      0: w[k] = sum A(i,k)*i
 *      1: w[k] = sum A(:,k)
 *      2: w[k] = mean (A(:,k))
 *      3: w[k] = min A(:,k)
 *      4: w[k] = max A(:,k)
 *      5: w[k] = A(k,k), if A(k,k) = 0 sum.
 * 
 * 
 * On success, this function returns a double* array containing the weights 
 */

double* SPEX_WAMF_get_weights
(
    SPEX_matrix* A,
    int64_t option
);


/* Purpose: A symbolic matrix addition. This gives a rough upper bound on 
 * the number of bits in C = A + B. 
 * 
 * On success, the matrix C is returned. On failure, NULL is returned
 * 
 */

SPEX_matrix* SPEX_WAMF_symbolic_add
(
    SPEX_matrix *A,      // Left matrix
    SPEX_matrix *B       // Right matrix
);


/* Purpose: x = x + A(:,j), where x is a dense vector and A(:,j) is sparse 
 * This function gives a rough upper bound of the number of bits present in x
 */

int64_t SPEX_WAMF_scatter
(
    SPEX_matrix *A, // Input matrix
    int64_t j,          // Column of A
    double beta,    // scale
    int64_t *w,         // workspace vector
    double *x,      // x = x + beta A(:,j)
    int64_t mark,       // location of C
    SPEX_matrix *C, // Set in C
    int64_t nz          // Number of nonzeros
);

/* Purpose: C = A*B. This function gives a rough upper bound in 
 * the number of bits in A*B. Note that A and B are not explicitly multiplied since
 * we are dealing with symbolic bit operations
 * 
 * On success, the matrix C is returned.
 * 
 */

SPEX_matrix* SPEX_WAMF_symbolic_multiply
(
    SPEX_matrix *A,  // Left matrix
    SPEX_matrix *B   // Right matrix
);


/* Purpose: drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 
 */

int64_t SPEX_WAMF_fkeep 
(
    SPEX_matrix *A, 
    int64_t (*fkeep) (int64_t, int64_t, double, void *), 
    void *other
);

/* Purpose: depth-first search and postorder of a tree rooted at node j */
int64_t SPEX_WAMF_tdfs 
(
    int64_t j, 
    int64_t k, 
    int64_t *head, 
    int64_t *next, 
    int64_t *post, 
    int64_t *stack
);

/* Purpose: clear the workspace array w in minimum degree ordering */
int64_t SPEX_WAMF_wclear 
(
    int64_t mark, 
    int64_t lemax, 
    int64_t *w, 
    int64_t n
);

/* Purpose: drop diagonal entries in a matrix when used as a function */

int64_t SPEX_WAMF_diag
(
    int64_t i, 
    int64_t j, 
    double aij, 
    void *other
);


/* Purpose: Expand a symbolic matrix*/
SPEX_info SPEX_WAMF_sparse_realloc
(
    SPEX_matrix* A
);

/* Purpose: Perform the wamf ordering */
int64_t* SPEX_WAMF_wamf
(
    int64_t order,   
    SPEX_matrix *A,   
    int64_t option,
    double alpha,
    bool aggressive
);

#endif
