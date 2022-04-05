//------------------------------------------------------------------------------
// SPEX_Update/SPEX_Update_Solve_mod: find the exact solution for Ax=b
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
 * x_handle: On input, it is a dense matrix containing the set of RHS vectors.
 *           On output, it is overwriten as the solution to LD^(-1)U x_out =
 *           x_in.
 *
 * F:        SPEX LU factorization in the updatable format
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

SPEX_info SPEX_Update_LU_Solve_mod // solves LD^(-1)U x_out = x_in
(
    SPEX_matrix **x_handle, // Input as a n*m dense matrix contains the
                            // right-hand-side vectora and output as the
                            // solution to the system. 
    const SPEX_factorization *F,// The SPEX LU factorization in dynamic_CSC
                            // format.
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
