//------------------------------------------------------------------------------
// SPEX_Update/SPEX_Update_Solve: find the exact solution for Ax=b
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function solves the linear system L(P,:)D^(-1)U(:,Q) x = b. It
 * essnetially serves as a wrapper for all forward and backward substitution
 * routines. This function always returns the solution matrix x as a mpz_t
 * matrix with additional pending scale factor. If a user desires to have
 * each value for each entries with pending scale factor applied, simply
 * compute x->v[j]->x[i]/x->scale and convert to double or mpfr output as
 * desired.
 *
 * Input/output arguments:
 *
 * x_handle: A pointer to the solution vectors. If x_handle = &b, then
 *           solution of the system will overwrite original b and return it
 *           as *x_handle. Otherwise, memory space will be allocated for
 *           x_handle to store the exact solution of the system
 *
 * b:        Set of RHS vectors
 *
 * L:        Lower triangular matrix. Unmodified on input/output
 *
 * U:        Upper triangular matrix. Unmodified on input/output
 *
 * A_scale:  Scale of the input matrix. Unmodified on input/output
 *
 * h:        array of n scalars, each entry should be >= -1
 *
 * rhos:     n-by-1 matrix that contains pivots
 *
 * P & Q_inv: the permutation vectors. unmodified on input/output.
 *
 * option:   command options
 */

#define SPEX_FREE_WORK                  \
    SPEX_FREE(v);                       \
    spex_delete_mpz_array(&x_col, n) ;

#define SPEX_FREE_ALL                   \
    SPEX_FREE_WORK                      \
    SPEX_matrix_free (&x, option) ;

#include "spex_update_internal.h"

SPEX_info SPEX_Update_Solve // solves Ax = b via REF LU factorization of A
(
    // Output
    SPEX_matrix **x_handle,    // solution to the system
    // input:
    SPEX_matrix *b,            // right hand side vector
    const SPEX_matrix *L,      // lower triangular matrix
    const SPEX_matrix *U,      // upper triangular matrix
    const mpq_t A_scale,    // scale of the input matrix
    int64_t *h,             // history vector
    const SPEX_matrix *rhos,// array of scaled pivots
    const int64_t *P,       // row permutation
    const int64_t *Q_inv,   // inverse of column permutation
    const SPEX_options* option // Command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_info info ;
    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    SPEX_REQUIRE(L, SPEX_DYNAMIC_CSC, SPEX_MPZ);
    SPEX_REQUIRE(U, SPEX_DYNAMIC_CSC, SPEX_MPZ);
    SPEX_REQUIRE(rhos, SPEX_DENSE, SPEX_MPZ);
    SPEX_REQUIRE(b,    SPEX_DENSE, SPEX_MPZ);
    if (!x_handle || !h || !P || !Q_inv || rhos->m != L->m ||
        L->m != L->n || L->n != U->m || U->n != U->m || L->m != b->m)
    {
        return SPEX_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // Declare and initialize workspace
    //--------------------------------------------------------------------------

    int64_t i, j, n = L->n;
    SPEX_matrix *x = NULL;   // final solution
    mpz_t *x_col = NULL;     // used to permute each col of x
    SPEX_vector *v = NULL;

    // allocate space for x_col and v
    x_col = spex_create_mpz_array(n);
    v = (SPEX_vector*) SPEX_malloc(sizeof(SPEX_vector));// v->x will be shallow
    if (!x_col || !v)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }
    v->x = NULL;
    v->i = NULL;
    v->nzmax = n;
    v->nz = n;
    SPEX_MPQ_SET_NULL(v->scale);// will not be used

    // allocate space for x if needed
    bool keep_b = (*x_handle != b);      // indicate if b will be reused
    if (keep_b)
    {
        *x_handle = NULL;
        // make a copy of b as x
        SPEX_CHECK(SPEX_matrix_copy(&x, SPEX_DENSE, SPEX_MPZ, b, option));
    }
    else
    {
        x = b;
        b = NULL;
    }

    for (j = 0; j < x->n; j++)
    {
        // TODO test with 2 columns b
        //----------------------------------------------------------------------
        // solve each column of b seperately
        //----------------------------------------------------------------------
        v->x = &(x->x.mpz[n*j]);

        // solve y for LD^(-1)y(P)=b, via forward substitution
        SPEX_CHECK(spex_update_forward_sub(v, L, P, rhos, h));

        // solve x for Ux(Q_inv) = y(P), via backward substitution
        SPEX_CHECK(spex_update_backward_sub(v, U, rhos, P, Q_inv));

        // permute x using P and Q_inv
        for (i = 0; i < n; i++)
        {
            SPEX_CHECK(SPEX_mpz_swap(x_col[i], v->x[P[Q_inv[i]]]));
        }
        // swap x->v[j]->x and x_col, then x->v[j]->x is permuted
        // and x_col unchanged
        for (i = 0; i < n; i++)
        {
            SPEX_CHECK(SPEX_mpz_swap(x_col[i], v->x[i]));
        }
    }

    //--------------------------------------------------------------------------
    // update the scale for the solution.
    //--------------------------------------------------------------------------
    // set the scaling factor x->scale *= rhos[n-1] / A_scale
    // the real solution is obtained by x->v[j]->x[i]/x->scale
    SPEX_CHECK(SPEX_mpz_mul(SPEX_MPQ_NUM(x->scale),
                            SPEX_MPQ_NUM(x->scale), rhos->x.mpz[n-1]));
    SPEX_CHECK(SPEX_mpq_canonicalize(x->scale));
    SPEX_CHECK(SPEX_mpq_div(x->scale, x->scale, A_scale));

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SPEX_FREE_WORK ;
    (*x_handle) = x ;
    return (SPEX_OK) ;
}
