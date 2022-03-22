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
 * UT:       Transpose of upper triangular matrix. Unmodified on input/output
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
    SPEX_vector_free(&v, option);

#define SPEX_FREE_ALL                   \
    SPEX_FREE_WORK                      \
    if (overwrite_b)  {(*x_handle) = x;}\
    else {SPEX_matrix_free (&x, option);}

#include "spex_update_internal.h"

SPEX_info spex_update_solve_internal// solves Ax = b via REF LU factorization of A
(
    // Output
    SPEX_matrix **x_handle, // a n*m dense matrix contains the solution to
                            // the system. If users wish to overwrite the
                            // solution to the right-hand-side matrix b, this
                            // can be provided as &b. Otherwise, new space will
                            // be allocated for x_handle
    // input:
    SPEX_matrix *b,         // a n*m dense matrix contains the right hand
                            // side vector
    const SPEX_matrix *L,   // a n*n dynamic_CSC matrix that gives the lower
                            // triangular matrix
    const SPEX_matrix *UT,  // a n*n dynamic_CSC matrix that gives the transpose
                            // of the upper triangular matrix
    const mpq_t A_scale,    // scale of the input matrix
    int64_t *h,             // history vector// TODO create a wrapper without h?
    const SPEX_matrix *rhos,// a n*1 dense matrix that gives the array of pivots
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

    SPEX_REQUIRE(L,  SPEX_DYNAMIC_CSC, SPEX_MPZ);
    SPEX_REQUIRE(UT, SPEX_DYNAMIC_CSC, SPEX_MPZ);
    SPEX_REQUIRE(rhos, SPEX_DENSE, SPEX_MPZ);
    SPEX_REQUIRE(b,    SPEX_DENSE, SPEX_MPZ);
    if (!x_handle || !h || !P || !Q_inv ||
        rhos->m != L->m || L->n  != L->m ||
        UT->m   != L->m || UT->n != L->m || L->m != b->m)
    {
        return SPEX_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // Declare and initialize workspace
    //--------------------------------------------------------------------------

    int64_t i, j, n = L->n;
    bool overwrite_b = (*x_handle == b); // indicate if b will be overwriten
    SPEX_matrix *x = NULL;               // final solution
    SPEX_vector *v = NULL;

    // allocate space for v and initialize
    v = (SPEX_vector*) SPEX_malloc(sizeof(SPEX_vector));
    if (!v)   { return SPEX_OUT_OF_MEMORY; }
    v->x = NULL;
    v->i = NULL;
    v->nzmax = n;
    v->nz = n;
    SPEX_MPQ_SET_NULL(v->scale);         // will not be used

    // allocate space for v->x
    v->x = spex_create_mpz_array(n);
    if (v->x == NULL)
    {
        SPEX_FREE(v);
        return SPEX_OUT_OF_MEMORY;
    }

    //--------------------------------------------------------------------------
    // allocate space for x if needed
    //--------------------------------------------------------------------------
    if (overwrite_b)
    {
        x = b;
        b = NULL;
    }
    else
    {
        *x_handle = NULL;
        // make a copy of b as x
        SPEX_CHECK(SPEX_matrix_copy(&x, SPEX_DENSE, SPEX_MPZ, b, option));
    }

    //--------------------------------------------------------------------------
    // solve each column of b seperately
    //--------------------------------------------------------------------------
    for (j = 0; j < x->n; j++)
    {
        // swap entries in j-th column of x with v
        for (i = 0; i < n; i++)
        {
            SPEX_CHECK(SPEX_mpz_swap(v->x[i], x->x.mpz[i+n*j]));
        }

        // solve y for LD^(-1)y(P)=b, via forward substitution
        SPEX_CHECK(spex_update_forward_sub(v, L, P, rhos, h));

        // solve x for Ux(Q_inv) = y(P), via backward substitution
        SPEX_CHECK(spex_update_backward_sub(v, UT, rhos, P, Q_inv));

        // permute v->x using P and Q_inv, and swap with x
        for (i = 0; i < n; i++)
        {
            SPEX_CHECK(SPEX_mpz_swap(x->x.mpz[i+n*j], v->x[P[Q_inv[i]]]));
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
