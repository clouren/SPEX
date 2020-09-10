//------------------------------------------------------------------------------
// SPEX_Chol/spex_chol_internal: include file for internal use in SPEX_Cholesky
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
// Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

// This file is not intended to be #include'd in user applications.  Use
// SPEX_Chol.h instead.

#ifndef SPEX_CHOL_INTERNAL_H
#define SPEX_CHOL_INTERNAL_H

#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wunused-value"

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------C Libraries------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Standard C libraries
#include <setjmp.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <inttypes.h>

// SuiteSparse headers
#include "SuiteSparse_config.h"
#include "colamd.h"
#include "amd.h"

//------------------------------------------------------------------------------
// debugging
//------------------------------------------------------------------------------

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

//------------------------------------------------------------------------------
//-------------------------Common Macros----------------------------------------
//------------------------------------------------------------------------------

#ifdef MATLAB_MEX_FILE
#include "mex.h"
#include "matrix.h"
#endif

#include "SPEX_Chol.h"

// Various other macros are inherited from SPEX_Util.h and may be seen in that file


// Free a pointer and set it to NULL.
#define SPEX_FREE(p)                        \
{                                           \
    SPEX_free (p) ;                         \
    (p) = NULL ;                            \
}                                           \

#ifndef FREE_WORKSPACE
#define FREE_WORKSPACE  
#endif

#define SPEX_MAX(a,b) (((a) > (b)) ? (a) : (b))

#define SPEX_MIN(a,b) (((a) < (b)) ? (a) : (b))

#define SPEX_FLIP(i) (-(i)-2)

#define SPEX_UNFLIP(i) (((i) < 0) ? SPEX_FLIP(i) : (i))

#define SPEX_MARKED(Ap,j) (Ap [j] < 0)

#define SPEX_MARK(Ap,j) { Ap [j] = SPEX_FLIP (Ap [j]) ; }

// Size of mpz_t, mpq_t and mpfr_t values
#define SIZE_MPZ  sizeof(mpz_t)

#define SIZE_MPQ  sizeof(mpq_t)

#define SIZE_MPFR sizeof(mpfr_t)


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
 * prevent error when SPEX_MP*_CLEAR is called before the variable is
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


    
    
//------------------------------------------------------------------------------
// Macros to utilize the default if option is NULL
//------------------------------------------------------------------------------

#define SPEX_DEFAULT_TOL 1

// Check parameter. If this = 1 then the solution to the system is checked
// for accuracy
#define SPEX_DEFAULT_CHECK false

// Column ordering used.
//  SPEX_NO_ORDERING = 0,           None: Not recommended for sparse matrices
//  SPEX_COLAMD = 1,                COLAMD: Default
//  SPEX_AMD = 2                    AMD
#define SPEX_DEFAULT_ORDER SPEX_COLAMD

// Defines printing to be done
#define SPEX_DEFAULT_PRINT_LEVEL 0

// MPFR precision used (quad is default)
#define SPEX_DEFAULT_PRECISION 128
    
#define SPEX_OPTION(option,parameter,default_value) \
    ((option == NULL) ? (default_value) : (option->parameter))

#define SPEX_OPTION_TOL(option) \
    SPEX_OPTION (option, tol, SPEX_DEFAULT_TOL)

#define SPEX_OPTION_CHECK(option) \
    SPEX_OPTION (option, check, false)

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
    
static inline int compare (const void * a, const void * b)
{
    return ( *(int64_t*)a - *(int64_t*)b );
}

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

// SPEX_LU uses SuiteSparse_config.printf_func instead of a mere call to printf
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

// ============================================================================
//                           Internal Functions
// ============================================================================

/* Purpose: Compute the elimination tree of A */
SPEX_info spex_Chol_etree 
(
    int64_t** tree,
    SPEX_matrix* A // Input matrix (must be SPD)
);

/* Purpose: This function computes the reach of the kth row of A onto the graph of L using the 
   elimination tree. This is more efficient than the SPEX_reach function 
   It finds the nonzero pattern of row k of L and uses the upper triangular 
   part of A(:,k) */
SPEX_info spex_Chol_ereach 
(
    int64_t* top_handle,
    SPEX_matrix *A,    // Matrix to be analyzed
    int64_t k,          // Node to start at
    int64_t* parent,    // ELimination Tree
    int64_t* s,         // Contains the nonzero pattern in s[top..n-1]
    int64_t* w          // Workspace array
);

/* Purpose: Depth-first search and postorder of a tree rooted at node j */
int64_t spex_Chol_tdfs 
(
    int64_t j,      // Root node
    int64_t k,      
    int64_t* head,  // Head of list
    int64_t* next,  // Next node in the list
    int64_t* post,  // Post ordered tree
    int64_t* stack  // Stack of nodes
);

/* Purpose: post order a forest */
SPEX_info spex_Chol_post 
(
    int64_t** post_handle,
    int64_t* parent,    // Parent[j] is parent of node j in forest
    int64_t n           // Number of nodes in the forest
);


/* Purpose: consider A(i,j), node j in ith row subtree and return lca(jprev,j) 
   Used to determine Column counts of cholesky factor*/
SPEX_info spex_Chol_leaf 
(
    int64_t* lca_handle,
    int64_t i, 
    int64_t j, 
    int64_t* first, 
    int64_t* maxfirst, 
    int64_t* prevleaf,
    int64_t* ancestor, 
    int64_t* jleaf
);

/* Purpose: Something*/
SPEX_info spex_Chol_counts 
(
    int64_t** c_handle,
    SPEX_matrix *A, 
    int64_t *parent, 
    int64_t *post
);

/* Purpose: This function performs the symmetric sparse REF triangular solve. for uplooking
 * Cholesky factorization. i.e., 
 * (LD) x = A(1:k-1,k). 
 */
SPEX_info spex_Up_Chol_triangular_solve // performs the sparse REF triangular solve
(
    int64_t *top_output,        // Output the beginning of nonzero pattern
    SPEX_matrix* L,              // partial L matrix
    SPEX_matrix* A,              // input matrix
    int64_t k,                    // iteration of algorithm
    int64_t* xi,                  // nonzero pattern vector
    int64_t* parent,              // Elimination tree
    int64_t* c,                   // Column pointers
    SPEX_matrix* rhos,              // sequence of pivots
    int64_t* h,                   // history vector
    SPEX_matrix* x                  // solution of system ==> kth column of L and U
);


/* Purpose: This solves the system L'x = b for Cholesky factorization */
SPEX_info spex_Chol_ltsolve 
(
    SPEX_matrix *L,    // The lower triangular matrix
    SPEX_matrix *x      // Solution vector
);

/* Purpose: This function performs the SLIP Cholesky factorization. This factorization
 * is done via n iterations of the sparse REF triangular solve function. The
 * overall factorization is PAP = LDL
 */
SPEX_info spex_Chol_Pre_Left_Factor         // performs the Up looking Cholesky factorization
(
    SPEX_matrix* A,
    SPEX_matrix** L_handle,              // partial L matrix
    int64_t* xi,                  // nonzero pattern vector
    int64_t* parent,              // Elimination tree
    SPEX_Chol_analysis * S,           // stores guess on nnz and column permutation
    int64_t* c                   // Column pointers
);

/* Purpose: This function performs the symmetric sparse REF triangular solve. i.e., 
 * (LD) x = A(:,k). 
 */
SPEX_info spex_Left_Chol_triangular_solve // performs the sparse REF triangular solve
(
    int64_t *top_output,        // Output the beginning of nonzero pattern
    SPEX_matrix* L,              // partial L matrix
    SPEX_matrix* A,              // input matrix
    int64_t k,                    // iteration of algorithm
    int64_t* xi,                  // nonzero pattern vector
    SPEX_matrix* rhos,              // sequence of pivots
    int64_t* h,                   // history vector
    SPEX_matrix* x,                  // solution of system ==> kth column of L and U
    int64_t* parent,
    int64_t* c
);

SPEX_info spex_Chol_forward_sub
(
    SPEX_matrix *L,   // lower triangular matrix
    SPEX_matrix *x,        // right hand side matrix of size n*numRHS
    SPEX_matrix *rhos      // sequence of pivots used in factorization
);

#endif

