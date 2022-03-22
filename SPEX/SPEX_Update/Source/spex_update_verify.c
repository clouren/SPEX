//------------------------------------------------------------------------------
// SPEX_Update/spex_update_verify.c: verify if A=LD^(-1)U
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to verify if A=L(P,:)D^(-1)U(:,Q)

#define SPEX_FREE_ALL                 \
    SPEX_matrix_free(&b, option);     \
    SPEX_matrix_free(&x, option);     \
    SPEX_matrix_free(&b2, option);    \
    SPEX_MPQ_CLEAR(x_scale);

#include "spex_update_internal.h"

SPEX_info spex_update_verify
(
    const SPEX_factorization *F,// LU factorization of A
    const SPEX_matrix *A,     // Input matrix
    int64_t *h,            // history vector
    const SPEX_options *option// command options
)
{
    SPEX_info info;
    int64_t tmp, i, n = F->L->n;
    int sgn;
    mpq_t x_scale; SPEX_MPQ_SET_NULL(x_scale);
    SPEX_matrix *b = NULL; // the dense right-hand-side matrix to be generated
    SPEX_matrix *x = NULL; // the dense solution matrix to be generated
    SPEX_matrix *b2 = NULL; // the dense matrix to store the result of A*x

    SPEX_CHECK(SPEX_mpq_init(x_scale));
    SPEX_CHECK(SPEX_mpq_set_ui(x_scale, 1, 1));
    SPEX_CHECK(SPEX_matrix_allocate(&b , SPEX_DENSE, SPEX_MPZ, n, 1, n, false,
        true, option));
    SPEX_CHECK(SPEX_matrix_allocate(&b2, SPEX_DENSE, SPEX_MPZ, n, 1, n, false,
        true, option));

    // -------------------------------------------------------------------------
    // generate random right-hand-size vector
    // -------------------------------------------------------------------------
    // initialize random number generator
    int seed = 10;
    srand(seed);
    for (i = 0; i < n; i++)
    {
        tmp = i+1;//rand(); //TODO
        SPEX_CHECK(SPEX_mpz_set_si(b->x.mpz[i], tmp));
    }

    // -------------------------------------------------------------------------
    // solve LD^(-1)Ux = b for x
    // -------------------------------------------------------------------------
    SPEX_CHECK(SPEX_Update_LU_Solve(&x, b, F, option));

    // -------------------------------------------------------------------------
    // compute b2 = A*x
    // -------------------------------------------------------------------------
    for (i = 0; i < n; i++)
    {
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, x->x.mpz[i]));
        if (sgn == 0) { continue;}

        for (int64_t p = 0; p < A->v[i]->nz; p++)
        {
            int64_t j = A->v[i]->i[p];
            // b2[j] += x[i]*A(j,i)
            SPEX_CHECK(SPEX_mpz_addmul(b2->x.mpz[j],
                x->x.mpz[i], A->v[i]->x[p]));
        }
    }
    // update b2->scale = x->scale*A->scale
    SPEX_CHECK(SPEX_mpq_mul(b2->scale, x->scale, A->scale));

    // -------------------------------------------------------------------------
    // check if b2 == b
    // -------------------------------------------------------------------------
    // set b2->scale = b2->scale/b->scale since we only want to compare the
    // integer values b2*b->scale and b*b->scale. b*b->scale are the values
    // stored in b->v->x. It can be shown that the resulted
    // b2->scale = rhos[n-1]
    SPEX_CHECK(SPEX_mpq_div(b2->scale, b2->scale, b->scale));
#ifdef SPEX_DEBUG
    SPEX_CHECK(SPEX_mpq_cmp_z(&sgn, b2->scale, F->rhos->x.mpz[n-1]));
    ASSERT(sgn == 0);
#endif
    for (i = 0; i < n; i++)
    {
        // set b2[i] = b2[i]/rhos[n-1] = b2[i]/b2->scale
        // This division will be exact since b are integers
        SPEX_CHECK(SPEX_mpz_divexact(b2->x.mpz[i],
                                     b2->x.mpz[i], SPEX_MPQ_NUM(b2->scale)));

        SPEX_CHECK(SPEX_mpz_cmp(&sgn, b2->x.mpz[i], b->x.mpz[i]));
        if (sgn != 0)
        {
            info = SPEX_INCORRECT;
            break;
        }
    }

    //--------------------------------------------------------------------------
    // Print info
    //--------------------------------------------------------------------------

    int pr = SPEX_OPTION_PRINT_LEVEL (option) ;
    if (info == SPEX_OK)
    {
        SPEX_PR1 ("Factorization is verified to be correct and exact.\n") ;
    }
    else if (info == SPEX_INCORRECT)
    {
        // This can never happen.
        SPEX_PR1 ("ERROR! Factorization is wrong. This is a bug; please "
                  "contact the authors of SPEX.\n") ;
    }

    SPEX_FREE_ALL;
    return info;
}
