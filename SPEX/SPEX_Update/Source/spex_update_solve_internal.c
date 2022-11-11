//------------------------------------------------------------------------------
// SPEX_Update/spex_update_solve_internal: find the exact solution for Ax=b
// with the the updatable LU factorization of A.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/* Purpose: This function is called by SPEX_update_solve or SPEX_update_tsolve
 * to solve the linear system Ax = b or ATx=b.
 *
 * Input/output arguments:
 *
 * x_handle: A pointer to the solution vectors. Memory space will be allocated
 *           for x_handle to store the exact solution of the system
 *
 * F:        SPEX LU or Cholesky factorization
 *
 * b:        Set of RHS vectors
 *
 * option:   command options
 */

#define SPEX_FREE_WORK                  \
{                                       \
    SPEX_FREE(h);                       \
    SPEX_vector_free(&v, option);       \
}

#define SPEX_FREE_ALL                   \
{                                       \
    SPEX_FREE_WORK                      \
    SPEX_matrix_free (&x, option);      \
}

#include "spex_update_internal.h"

SPEX_info spex_update_solve_internal
(
    // Output
    SPEX_matrix *x_handle, // a m*n dense matrix contains the solution to
                            // the system.
    // input:
    SPEX_factorization *F,  // The SPEX LU or Cholesky factorization
    const SPEX_matrix b,   // a m*n dense matrix contains the right-hand-side
                            // vector
    const bool transpose,   // whether computing Ax=b or ATx=b
    const SPEX_options* option // Command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_info info ;
    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    SPEX_REQUIRE(b,    SPEX_DENSE, SPEX_MPZ);
    if (!x_handle || !F || (F->kind != SPEX_LU_FACTORIZATION &&
        F->kind != SPEX_CHOLESKY_FACTORIZATION) || F->L->m != b->m)
    {
        return SPEX_INCORRECT_INPUT;
    }
    *x_handle = NULL;
    
    // convert F to updatable format
    info = SPEX_factorization_convert(F, true, option);
    if (info != SPEX_OK) return info;

    //--------------------------------------------------------------------------
    // Declare and initialize workspace
    //--------------------------------------------------------------------------

    SPEX_matrix L, UT, rhos = F->rhos;
    int64_t *P, *Q_inv;
    int64_t i, j, n = F->L->n;
    int64_t *h = NULL;                   // history vector
    SPEX_vector *v = NULL;               // temp mpz vector 
    SPEX_matrix x = NULL;               // final solution

    if (F->kind == SPEX_CHOLESKY_FACTORIZATION)
    {
        L = F->L;
        UT = F->L;
        P = F->P_perm;
        Q_inv = F->Pinv_perm;
    }
    else // F->kind == SPEX_LU_FACTORIZATION
    {
        if (!transpose)
        {
            L = F->L;
            UT = F->U;
            P = F->P_perm;
            Q_inv = F->Qinv_perm;
        }
        else
        {
            L = F->U;
            UT = F->L;
            P = F->Q_perm;
            Q_inv = F->Pinv_perm;
        }
    }

    // allocate space for v and initialize
    v = (SPEX_vector*) SPEX_malloc(sizeof(SPEX_vector));
    if (!v)   { return SPEX_OUT_OF_MEMORY; }
    v->x = NULL;
    v->i = NULL; // leave this NULL since v is dense
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
    SPEX_CHECK(SPEX_matrix_allocate(&x, SPEX_DENSE, SPEX_MPQ, b->m, b->n, 0,
        false, true, option));

    // compute the scale first with scale = b->scale * rhos[n-1] / A_scale
    // the real solution is obtained by x->v[j]->x[i]/x->scale
    SPEX_CHECK(SPEX_mpq_set_z(x->scale, rhos->x.mpz[n-1]));
    SPEX_CHECK(SPEX_mpq_mul(x->scale, x->scale, b->scale));
    SPEX_CHECK(SPEX_mpq_div(x->scale, x->scale, F->scale_for_A));

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

        // permute v->x using P and Q_inv, and apply x->scale
        for (i = 0; i < n; i++)
        {
            int64_t p = i+n*j;
            SPEX_CHECK(SPEX_mpq_set_z(x->x.mpq[p], v->x[P[Q_inv[i]]]));
            SPEX_CHECK(SPEX_mpq_div(x->x.mpq[p], x->x.mpq[p], x->scale));
        }
    }

    // reset the scale for the solution.
    SPEX_CHECK(SPEX_mpq_set_ui(x->scale, 1, 1));

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SPEX_FREE_WORK ;
    (*x_handle) = x ;
    return (SPEX_OK) ;
}
