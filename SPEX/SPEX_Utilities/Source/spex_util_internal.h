//------------------------------------------------------------------------------
// SPEX_Utilities/spex_util_internal: include file for internal use in SPEX_Utility functions
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// This file is not intended to be #include'd in user applications.  Use
// SPEX.h instead.

#ifndef SPEX_UTIL_INTERNAL_H
#define SPEX_UTIL_INTERNAL_H

#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-value"

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------C Libraries------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Standard C libraries
#include <setjmp.h>
#include <stdarg.h>

// SuiteSparse headers
#include "SuiteSparse_config.h"

// macros for MPZ/MPQ/MPFR struct
#include "spex_gmp.h"

// user-callable functions
#include "SPEX.h"


//------------------------------------------------------------------------------
// debugging
//------------------------------------------------------------------------------

#undef SPEX_DEBUG
// uncomment this line to enable debugging
 #define SPEX_DEBUG

#ifdef SPEX_DEBUG

#ifdef MATLAB_MEX_FILE

#define ASSERT(x)                                                             \
{                                                                             \
    if (!(x))                                                                 \
    {                                                                         \
        mexErrMsgTxt ("assertion failed: %s line %d\n", __FILE__, __LINE__) ; \
    }                                                                         \
}

#else

#include <assert.h>
#define ASSERT(x) assert (x)

#endif

#else

// debugging disabled
#define ASSERT(x)

#endif

#undef  GOTCHA
#define GOTCHA                                                              \
{                                                                           \
    printf ("GOTCHA: %s, line %d\n", __FILE__, __LINE__) ;                  \
    abort ( ) ;                                                             \
}

//------------------------------------------------------------------------------
//-------------------------Common Macros----------------------------------------
//------------------------------------------------------------------------------

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#endif

#define SPEX_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define SPEX_MIN(a,b) (((a) < (b)) ? (a) : (b))
#define SPEX_FLIP(i) (-(i)-2)
#define SPEX_UNFLIP(i) (((i) < 0) ? SPEX_FLIP(i) : (i))
#define SPEX_MARKED(Ap,j) (Ap [j] < 0)
#define SPEX_MARK(Ap,j) { Ap [j] = SPEX_FLIP (Ap [j]) ; }

// SPEX_CHECK(method) is a macro that calls a SPEX method and checks the
// status; if a failure occurs, it frees all allocated workspace and returns
// the error status to the caller.  To use SPEX_CHECK, the #include'ing file
// must declare a SPEX_info info, and must define SPEX_FREE_ALL as a macro that
// frees all workspace if an error occurs. The method can be a scalar as well,
// so that SPEX_CHECK(info) works.

// To create a SPEX_FREE_WORKSPACE and SPEX_FREE_ALL, use this approach:

/*
    #define SPEX_FREE_WORKSPACE \
    {                           \
        free all workspace      \
    }
    #define SPEX_FREE_ALL       \
    {                           \
        SPEX_FREE_WORKSPACE ;   \
        free all output objects \
    }
    #include "this file.h"
*/

// then, since SPEC_FREE_ALL is already defined, it isn't #defined below.
//
// Definitions of these macros:
//
//  SPEX_FREE_ALL:  frees everything, including both workspace and the result
//  SPEX_FREE_WORKSPACE: frees just the internal workspace of the method, not
//  the final result

// SPEX_FREE_ALL: Frees workspace and memory that would be returned. Only used
// when there is an error and we need to cleanly close things (like out of
// memory conditions and failed SPEX_CHECK)

#ifndef SPEX_FREE_ALL
#define SPEX_FREE_ALL
#endif

// SPEX_FREE_WORKSPACE: Frees all workspace (but not memory that was allocated
// and meant to be returned). Used usually at the end of the function.  This
// workspace mechanism should be used even if the function only has one
// workspace variable (it can be risky otherwise with further development)

#ifndef SPEX_FREE_WORKSPACE
#define SPEX_FREE_WORKSPACE
#endif

// Local variables (only declared, allocated and freed inside an if, for
// example) do not go inside the workspace.

#ifdef SPEX_DEBUG

    #define SPEX_CHECK(method)      \
    {                               \
        info = (method) ;           \
        if (info != SPEX_OK)        \
        {                           \
            printf("file %s line %d\n",__FILE__,__LINE__);\
            SPEX_FREE_ALL ;         \
            return (info) ;         \
        }                           \
    }

#else

    #define SPEX_CHECK(method)      \
    {                               \
        info = (method) ;           \
        if (info != SPEX_OK)        \
        {                           \
            SPEX_FREE_ALL ;         \
            return (info) ;         \
        }                           \
    }

#endif

//------------------------------------------------------------------------------
// printing control
//------------------------------------------------------------------------------

// SPEX uses SuiteSparse_config.printf_func instead of a mere call to printf
// (the default function is printf, or mexPrintf when in MATLAB).  If this
// function pointer is NULL, no printing is done.

#define SPEX_PRINTF(...)                                    \
{                                                           \
    if (SuiteSparse_config.printf_func != NULL)             \
    {                                                       \
        SuiteSparse_config.printf_func (__VA_ARGS__) ;      \
    }                                                       \
}

#define SPEX_PR1(...) { if (pr >= 1) SPEX_PRINTF (__VA_ARGS__) }
#define SPEX_PR2(...) { if (pr >= 2) SPEX_PRINTF (__VA_ARGS__) }
#define SPEX_PR3(...) { if (pr >= 3) SPEX_PRINTF (__VA_ARGS__) }



//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-----------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// TODO: What should we do to fix this?

// Tolerance used in the pivoting schemes. This number can be anything in
// between 0 and 1. A value of 0 selects the diagonal element exclusively and a
// value of 1 selects the smallest or largest pivot exclusively only in a
// tolerance pivoting based method
#define SPEX_DEFAULT_TOL 1


// Pivoting scheme used for SPEX Left LU.
//  SPEX_SMALLEST = 0,              Smallest pivot
//  SPEX_DIAGONAL = 1,              Diagonal pivoting
//  SPEX_FIRST_NONZERO = 2,         First nonzero per column chosen as pivot
//  SPEX_TOL_SMALLEST = 3,          Diagonal pivoting with tolerance for small
//  SPEX_TOL_LARGEST = 4,           Diagonal pivoting with tolerance for large
//  SPEX_LARGEST = 5                Largest pivot
#define SPEX_DEFAULT_PIVOT SPEX_SMALLEST

// Column ordering used.
//  SPEX_NO_ORDERING = 0,           None: Not recommended for sparse matrices
//  SPEX_COLAMD = 1,                COLAMD: Default
//  SPEX_AMD = 2                    AMD
#define SPEX_DEFAULT_ORDER SPEX_DEFAULT_ORDERING

// Defines printing to be done
#define SPEX_DEFAULT_PRINT_LEVEL 0

// MPFR precision used (quad is default)
#define SPEX_DEFAULT_PRECISION 128

// Defines the algorithm used
// SPEX_ALGORITHM_DEFAULT = 0,    Defaults: Left for LU, Up for Chol, Gram for QR looking LU factorization //TODO maybe rename default...
// SPEX_LU_LEFT = 1,              Left looking LU factorization
// SPEX_CHOL_LEFT = 2,            Left looking Cholesky factorization
// SPEX_CHOL_UP = 3,              Up looking Cholesky factorization
// SPEX_QR_GRAM = 4               Default factorization for QR
#define SPEX_DEFAULT_ALGORITHM SPEX_DEFAULT

//------------------------------------------------------------------------------
// Type of MPFR rounding used.
//------------------------------------------------------------------------------

// The MPFR library utilizes an internal rounding scheme. The options are
//  MPFR_RNDN: round to nearest (roundTiesToEven in IEEE 754-2008),
//  MPFR_RNDZ: round toward zero (roundTowardZero in IEEE 754-2008),
//  MPFR_RNDU: round toward plus infinity (roundTowardPositive in
//             IEEE 754-2008),
//  MPFR_RNDD: round toward minus infinity (roundTowardNegative in
//             IEEE 754-2008),
//  MPFR_RNDA: round away from zero.
//  MPFR_RNDF: faithful rounding. This is not stable
//
// SPEX utilizes MPFR_RNDN by default.

#define SPEX_DEFAULT_MPFR_ROUND MPFR_RNDN

//------------------------------------------------------------------------------
// Macros to utilize the default if option is NULL
//------------------------------------------------------------------------------

#define SPEX_OPTION(option,parameter,default_value) \
    ((option == NULL) ? (default_value) : (option->parameter))

#define SPEX_OPTION_TOL(option) \
    SPEX_OPTION (option, tol, SPEX_DEFAULT_TOL)

#define SPEX_OPTION_PIVOT(option) \
    SPEX_OPTION (option, pivot, SPEX_DEFAULT_PIVOT)

#define SPEX_OPTION_ORDER(option) \
    SPEX_OPTION (option, order, SPEX_DEFAULT_ORDER)

#define SPEX_OPTION_PREC(option) \
    SPEX_OPTION (option, prec, SPEX_DEFAULT_PRECISION)

#define SPEX_OPTION_PRINT_LEVEL(option) \
    SPEX_OPTION (option, print_level, SPEX_DEFAULT_PRINT_LEVEL)

#define SPEX_OPTION_ROUND(option) \
    SPEX_OPTION (option, round, SPEX_DEFAULT_MPFR_ROUND)

#define SPEX_OPTION_ALGORITHM(option) \
    SPEX_OPTION (option, algo, SPEX_DEFAULT_ALGORITHM)



// ============================================================================
//                           Internal Functions
// ============================================================================

// check if SPEX_initialize* has been called
bool spex_initialized ( void ) ;        // true if called, false if not
void spex_set_initialized (bool s) ;    // set global initialzed flag to s


//------------------------------------------------------------------------------
// spex_create_mpfr_array: create a 1D mpfr_t array
//------------------------------------------------------------------------------

// Creates a simple 1D array, where A[i] is an entry of type mpfr_t.

/* Purpose: This function creates a MPFR array of desired precision*/
mpfr_t* spex_create_mpfr_array
(
    int64_t n,     // size of the array
    const SPEX_options* option
);

//------------------------------------------------------------------------------
// spex_create_mpq_array: Creates a 1D array, whose entries are all mpq_t type.
//------------------------------------------------------------------------------

// Creates a simple 1D array, where A[i] is an entry of type mpq_t.

/* Purpose: This function creates an mpq array of size n.
 * This function must be called for all mpq arrays created.
 */
mpq_t* spex_create_mpq_array
(
    int64_t n              // size of the array
);

SPEX_info spex_create_mpq (mpq_t x);
//------------------------------------------------------------------------------
// spex_create_mpz_array: create a 1D mpz_t array
//------------------------------------------------------------------------------

// Creates a simple 1D array, where A[i] is an entry of type mpz_t.

/* Purpose: This function creates an mpz array of size n and allocates
 * default size.
 */
mpz_t* spex_create_mpz_array
(
    int64_t n              // Size of x
);


//------------------------------------------------------------------------------
// spex_delete_mpz_array: delete a 1D mpz_t array
//------------------------------------------------------------------------------

// Delete a simple 1D array, where A[i] is an entry of type mpz_t.

/* Purpose: This function deletes a mpz array of size n
 */
void spex_delete_mpz_array
(
    mpz_t **x,      // mpz array to be deleted
    int64_t n       // Size of x
);



/* Purpose: This function converts a double array of size n to an appropriate
 * mpz array of size n. To do this, the number is multiplied by 10^17 then, the
 * GCD is found. This function allows the use of matrices in double precision
 * to work with SPEX.
 */
SPEX_info spex_expand_double_array
(
    mpz_t *x_out,   // integral final array
    double* x,      // double array that needs to be made integral
    mpq_t scale,    // the scaling factor used (x_out = scale * x)
    int64_t n,      // size of x
    const SPEX_options* option
);

/* Purpose: This function converts a mpfr array of size n and precision prec to
 * an appropriate mpz array of size n. To do this, the number is multiplied by
 * the appropriate power of 10 then the gcd is found. This function allows mpfr
 * arrays to be used within SPEX.
 */
SPEX_info spex_expand_mpfr_array
(
    mpz_t *x_out,   // integral final array
    mpfr_t* x,      // mpfr array to be expanded
    mpq_t scale,    // scaling factor used (x_out = scale*x)
    int64_t n,      // size of x
    const SPEX_options *option // command options containing the prec for mpfr
);

/* Purpose: This function converts a mpq array of size n into an appropriate mpz
 * array of size n. To do this, the lcm of the denominators is found as a
 * scaling factor. This function allows mpq arrays to be used in SPEX
 */
SPEX_info spex_expand_mpq_array
(
    mpz_t *x_out, // integral final array
    mpq_t* x,     // mpq array that needs to be converted
    mpq_t scale,  // scaling factor. x_out = scale*x
    int64_t n,     // size of x
    const SPEX_options* option // Command options
);

/* Purpose: This function converts a mpq matrix of size m*n into an appropriate
 * mpz matrix of size m*n. To do this, the lcm of the denominators is found as a
 * scaling factor. This function allows mpq matrix to be used in SPEX.
 * on output, x2 is modified.
 */
SPEX_info spex_expand_mpq_mat
(
    mpz_t **x_out,// integral final mat
    mpq_t **x,    // mpq mat that needs to be converted
    mpq_t scale,  // scaling factor. x_out = scale*x
    int64_t m,    // number of rows of x
    int64_t n     // number of columns of x
);


// typecast a double value to int64, accounting for Infs and Nans
static inline int64_t spex_cast_double_to_int64 (double x)
{
    if (isnan (x))
    {
        return (0) ;
    }
    else if (x > INT64_MAX)
    {
        return (INT64_MAX) ;
    }
    else if (x < INT64_MIN)
    {
        return (INT64_MIN) ;
    }
    else
    {
        return ((int64_t) (x)) ;
    }
}

SPEX_info spex_cast_array
(
    void *Y,                // output array, of size n
    SPEX_type ytype,        // type of Y
    void *X,                // input array, of size n
    SPEX_type xtype,        // type of X
    int64_t n,              // size of Y and X
    mpq_t y_scale,          // scale factor applied if y is mpz_t
    const mpq_t x_scale,          // scale factor applied if x is mpz_t
    const SPEX_options *option
) ;

SPEX_info spex_cast_matrix
(
    SPEX_matrix **Y_handle,     // nz-by-1 dense matrix to create
    SPEX_type Y_type,           // type of Y
    const SPEX_matrix *A,             // matrix with nz entries
    const SPEX_options *option
) ;

SPEX_info spex_CSC_mpz_to_dynamic
(
    SPEX_matrix **A_handle,       // converted SPEX_mat matrix
    const SPEX_matrix *B,         // original matrix
    const SPEX_options *option
);

SPEX_info spex_dynamic_to_CSC_mpz
(
    SPEX_matrix **A_handle,       // converted CSC matrix
    const SPEX_matrix *B,         // original matrix
    const int64_t nnz,            // number of nonzeros in B
    const SPEX_options *option
);


/* Purpose: This function collapses a SPEX matrix. Essentially it shrinks the
 * size of x and i. so that they only take up the number of elements in the
 * matrix. For example if A->nzmax = 1000 but nnz(A) = 500, i and x are of size
 * 1000, so this function shrinks them to size 500.
 */
SPEX_info spex_sparse_collapse
(
    SPEX_matrix* A // matrix to be shrunk
);

/* Purpose: This function expands a SPEX matrix by doubling its size. It
 * merely expands x and i and does not initialize/allocate the values.
 */
SPEX_info spex_sparse_realloc
(
    SPEX_matrix* A // the matrix to be expanded
);

/* Purpose: This function gets a copy of a row-wise permuted dense matrix as
 * A_out = P*A_in.
 */

SPEX_info spex_permute_dense_matrix
(
    SPEX_matrix **A_handle,     // permuted A
    const SPEX_matrix *A_in,    // unpermuted A (not modified)
    const int64_t *P,           // row permutation
    const SPEX_options* option
);

/* Purpose: This function multiplies matrix x a scalar
 */
SPEX_info spex_matrix_mul   // multiplies x by a scalar
(
    SPEX_matrix *x,         // matrix to be multiplied
    const mpz_t scalar      // scalar to multiply by
) ;

/* Purpose: p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1]
 * into c.  This function is lightly modified from CSparse.
 */
SPEX_info spex_cumsum
(
    int64_t *p,          // vector to store the sum of c
    int64_t *c,          // vector which is summed
    int64_t n            // size of c
);

/* Purpose: perform basic check for a given factorization
 */
SPEX_info spex_factorization_basic_check
(
    SPEX_factorization *F // The factorization to check
);

SPEX_info spex_colamd
(
    int64_t **perm_handle,
    int64_t *nnz,
    const SPEX_matrix *A,
    const SPEX_options* option
);

SPEX_info spex_amd
(
    int64_t **perm_handle,
    int64_t *nnz,
    const SPEX_matrix *A,
    const SPEX_options* option
);

// (void *) pointer to the values of A.  A must be non-NULL with a valid type
#define SPEX_X(A)                                                           \
    ((A->type == SPEX_MPZ  ) ? (void *) A->x.mpz   :                        \
    ((A->type == SPEX_MPQ  ) ? (void *) A->x.mpq   :                        \
    ((A->type == SPEX_MPFR ) ? (void *) A->x.mpfr  :                        \
    ((A->type == SPEX_INT64) ? (void *) A->x.int64 : (void *) A->x.fp64))))


// return an error if A->kind (csc, triplet, dense) is wrong
#define SPEX_REQUIRE_KIND(A,required_kind) \
    if (A == NULL || A->kind != required_kind) return (SPEX_INCORRECT_INPUT) ;

#define ASSERT_KIND(A,required_kind) \
    ASSERT (A != NULL && A->kind == required_kind)

// return an error if A->type (mpz, mpq, mpfr, int64, or double) is wrong
#define SPEX_REQUIRE_TYPE(A,required_type) \
    if (A == NULL || A->type != required_type) return (SPEX_INCORRECT_INPUT) ;

#define ASSERT_TYPE(A,required_type) \
    ASSERT (A != NULL && A->type == required_type)

// return an error if A->kind or A->type is wrong
#define SPEX_REQUIRE(A,required_kind,required_type)     \
    SPEX_REQUIRE_KIND (A,required_kind) ;               \
    SPEX_REQUIRE_TYPE (A,required_type) ;

#define ASSERT_MATRIX(A,required_kind,required_type)    \
    ASSERT_KIND (A,required_kind) ;                     \
    ASSERT_TYPE (A,required_type) ;

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

#endif
