//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_solve: find the exact solution for Ax=b
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/* Purpose: This function solves the linear system L(P,:)D^(-1)U(:,Q) x = b. It
 * essnetially serves as a wrapper for all forward and backward substitution
 * routines. This function always returns the solution matrix x as a mpz_t
 * matrix with additional pending scale factor. If a user desires to have
 * each value for each entries with pending scale factor applied, simply
 * compute x->v[j]->x[i]/x->scale and convert to double or mpfr output as
 * desired.
 *
 * Input/output arguments://TODO
 *
 * x_handle: A pointer to the solution vectors. Unitialized on input.
 *           on output, contains the exact solution of the system
 *
 * b:        Set of RHS vectors
 *
 * L:        Lower triangular matrix. Unmodified on input/output
 *
 * U:        Upper triangular matrix. Unmodified on input/output
 *
 * A_scale:  Scale of the input matrix. Unmodified on input/output
 *
 * S:        a dense 3*n matrix of scale factor but stored as a vector. S(1,:)
 *           is the pending scale for L, S(2,:) is the pending scale for U,
 *           and S(3, :) is the pending scale for both L and U.
 *
 * sd:       array of pivots with pending scale applied
 *
 * P & Q_inv: the permutation vectors. unmodified on input/output.
 *
 * keep_b:   indicate if caller wants to keep the vector b. If true, this
 *           function will make a copy of b before forward_sub and
 *           backward_sub. Otherwise, x->v[j]->x will be directly set as
 *           b->v[j]->x and b->v[j]->x  will be reset to NULL.
 *
 * option:   command options
 */

#define SPEX_FREE_WORK                  \
    SPEX_delete_mpz_array(&x_col, n) ;

#define SPEX_FREE_ALL                   \
    SPEX_FREE_WORK                      \
    SPEX_mat_free (&x) ;

#include "spex_lu_update_internal.h"

// TODO: fix to enable solving A'x = c
SPEX_info SPEX_solve        // solves Ax = b via REF LU factorization of A
(
    // Output
    SPEX_mat **x_handle,    // solution to the system
    // input:
    SPEX_mat *b,            // right hand side vector
    const SPEX_mat *L,      // lower triangular matrix
    const SPEX_mat *U,      // upper triangular matrix
    const mpq_t A_scale,    // scale of the input matrix
    const bool Is_trans,    // true if solving A'x = b
    int64_t *h,             // history vector
    SPEX_matrix *S,         // the pending scale factor matrix
    const mpz_t *sd,        // array of scaled pivots
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

    if (!x_handle || !b || !h || !L || !U || !S || !sd || !P || !Q_inv ||
        L->m != L->n || L->n != U->m || U->n != U->m)
    {
        return SPEX_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // Declare and initialize workspace
    //--------------------------------------------------------------------------

    int64_t i, j, n = L->n;
    SPEX_mat *x = NULL;      // final solution
    mpz_t *x_col = NULL;     // used to permute each col of x

    // allocate space for x_col
    x_col = SPEX_create_mpz_array(n);
    if (!x_col)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    // allocate space for x if needed
    bool keep_b = (*x_handle != b);      // indicate if b will be reused
    if (keep_b)
    {
        *x_handle = NULL;
        // even though x will be dense, we initialize it as sparse, which make
        // each of its columns initialized with length 0, and we will allocate
        // space for its mpz_t vector later depend on the value keep_b 
        SPEX_CHECK(SPEX_mat_alloc(&x, b->n, n, true));

        // x->scale = b->scale
        SPEX_CHECK(SPEX_mpq_set(x->scale, b->scale));
    }
    else
    {
        x = b;
        b = NULL;
    }


    for (j = 0; j < x->n; j++)
    {
        //----------------------------------------------------------------------
        // solve each column of b seperately
        //----------------------------------------------------------------------
        if (keep_b)
        {
            // remove x->v[j]->i which is not needed for dense vector
            SPEX_FREE(x->v[j]->i);
            // just need to allocate space for x, while let i still be NULL
            x->v[j]->x = spex_create_mpz_array(n);
            if (x->v[j]->x == NULL)
            {
                SPEX_FREE_WORK;
                return SPEX_OUT_OF_MEMORY;
            }

            // make a copy of b->v[j] to x->v[j]
            // Notice that we don't need to permute b here since rows of L are
            // maintained in the same order as those of A
            for (i = 0; i < n; i++)
            {
                SPEX_CHECK(SPEX_mpz_set(x->v[j]->x[i], b->v[j]->x[i]));
            }
            x->v[j]->nzmax = n;
        }

        // solve y for LD^(-1)y(P)=b, via forward substitution
        SPEX_CHECK(spex_forward_sub(x->v[j], L, P, sd, S, h, Is_trans));

        // solve x for Ux(Q_inv) = y(P), via backward substitution
        SPEX_CHECK(spex_backward_sub(x->v[j], U, (const SPEX_matrix*) S, sd,
            P, Q_inv, Is_trans));

        // permute x using P and Q_inv
        for (i = 0; i < n; i++)
        {
            SPEX_CHECK(SPEX_mpz_swap(x_col[i], x->v[j]->x[P[Q_inv[i]]]));
        }
        // swap x->v[j]->x and x_col, then x->v[j]->x is permuted
        // and x_col unchanged
        mpz_t *tmp = x->v[j]->x; x->v[j]->x = x_col; x_col = tmp;
    }

    //--------------------------------------------------------------------------
    // update the scale for the solution.
    //--------------------------------------------------------------------------
    // set the scaling factor x->scale *= sd[n-1] / A_scale
    // the real solution is obtained by x->v[j]->x[i]/x->scale
    SPEX_CHECK(SPEX_mpz_mul(SPEX_MPQ_NUM(x->scale),
                            SPEX_MPQ_NUM(x->scale), sd[n-1]));
    SPEX_CHECK(SPEX_mpq_canonicalize(x->scale));
    SPEX_CHECK(SPEX_mpq_div(x->scale, x->scale, A_scale));

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SPEX_FREE_WORK ;
    (*x_handle) = x ;
    return (SPEX_OK) ;
}
