//------------------------------------------------------------------------------
// SPEX_Update/SPEX_Update_LU_Solve: find the exact solution for Ax=b with the
// the updatable LU factorizaiton of A.
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
 * x_handle: A pointer to the solution vectors. Memory space will be allocated
 *           for x_handle to store the exact solution of the system
 *
 * b:        Set of RHS vectors
 *
 * F:        SPEX LU factorization in the updatable format
 *
 * option:   command options
 */

#define SPEX_FREE_WORK                  \
    SPEX_FREE(h);                       \
    SPEX_vector_free(&v, option);

#define SPEX_FREE_ALL                   \
    SPEX_FREE_WORK                      \
    SPEX_matrix_free (&x, option);

#include "spex_update_internal.h"

SPEX_info SPEX_Update_LU_Solve // solves Ax = b via REF LU factorization of A
(
    // Output
    SPEX_matrix **x_handle, // a m*n dense matrix contains the solution to
                            // the system.
    // input:
    const SPEX_matrix *b,         // a m*n dense matrix contains the right-hand-side
                            // vector
    const SPEX_factorization *F,// The SPEX LU factorization in updatable form
    const SPEX_options* option // Command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    SPEX_REQUIRE(b,    SPEX_DENSE, SPEX_MPZ);
    if (!x_handle || !F || F->L->m != b->m)
    {
        return SPEX_INCORRECT_INPUT;
    }
    *x_handle = NULL;

    //--------------------------------------------------------------------------
    // Declare and initialize workspace
    //--------------------------------------------------------------------------

    SPEX_info info ;
    SPEX_matrix *L = F->L, *UT = F->U, *rhos = F->rhos;
    int64_t *P = F->P_perm, *Q_inv = F->Qinv_perm;
    int64_t i, j, n = L->n;
    int64_t *h = NULL;                   // history vector
    SPEX_vector *v = NULL;               // temp mpz vector 
    SPEX_matrix *x = NULL;               // final solution

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

    // allocate space for h
    h = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    if (h == NULL)
    {
        SPEX_vector_free(&v, option);
        return SPEX_OUT_OF_MEMORY;
    }

    // make a copy of b as x
    //SPEX_CHECK(SPEX_matrix_copy(&x, SPEX_DENSE, SPEX_MPZ, b, option));
    // allocate space for x
    SPEX_CHECK(SPEX_matrix_allocate(&x, SPEX_DENSE, SPEX_MPZ, b->m, b->n, 0,
        false, true, option));
    SPEX_CHECK(SPEX_mpq_set(x->scale, b->scale));

    //--------------------------------------------------------------------------
    // solve each column of b seperately
    //--------------------------------------------------------------------------
    for (j = 0; j < x->n; j++)
    {
        // copy entries in j-th column of b with v
        for (i = 0; i < n; i++)
        {
            SPEX_CHECK(SPEX_mpz_set(v->x[i], b->x.mpz[i+n*j]));
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
    SPEX_CHECK(SPEX_mpq_div(x->scale, x->scale, F->scale_for_A));

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SPEX_FREE_WORK ;
    (*x_handle) = x ;
    return (SPEX_OK) ;
}
