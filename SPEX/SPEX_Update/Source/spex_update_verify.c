//------------------------------------------------------------------------------
// SPEX_Update/spex_update_verify: verify if A=LD^(-1)U
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2023, Christopher Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Purpose: This function is to verify if A=L(P,:)D^(-1)U(:,Q) after
// factorization update. This is done by solving LD^(-1)U*x=b via the updated
// factorization, and check if A*x=b holds rational-arthmetically. This
// function is only used for debugging purposes, as the routines
// within SPEX are gauranteed to be exact.

#define SPEX_FREE_ALL                 \
    SPEX_matrix_free(&b, option);     \
    SPEX_matrix_free(&x, option);     \
    SPEX_matrix_free(&b2, option);    \
    SPEX_MPQ_CLEAR(temp);

#include "spex_update_internal.h"

#ifdef SPEX_DEBUG

SPEX_info spex_update_verify
(
    bool *Is_correct,           // if factorization is correct
    SPEX_factorization F,       // LU factorization of A
    const SPEX_matrix A,        // Input matrix Dynamic_CSC MPZ
    const SPEX_options option   // command options
)
{

    SPEX_info info = SPEX_OK;
    int64_t tmp, i, n = F->L->n;
    int r;
    mpq_t temp; SPEX_MPQ_SET_NULL(temp);
    SPEX_matrix b = NULL; // the dense right-hand-side matrix to be generated
    SPEX_matrix x = NULL; // the dense solution matrix to be generated
    SPEX_matrix b2 = NULL; // the dense matrix to store the result of A*x

    SPEX_MPQ_INIT(temp);
    SPEX_CHECK(SPEX_matrix_allocate(&b , SPEX_DENSE, SPEX_MPZ, n, 1, n, false,
        true, option));
    SPEX_CHECK(SPEX_matrix_allocate(&b2, SPEX_DENSE, SPEX_MPQ, n, 1, n, false,
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
        SPEX_MPZ_SET_SI(b->x.mpz[i], tmp);
    }

    // -------------------------------------------------------------------------
    // solve LD^(-1)Ux = b for x
    // -------------------------------------------------------------------------
    SPEX_CHECK(SPEX_update_solve(&x, F, b, option));

    // -------------------------------------------------------------------------
    // compute b2 = A*x
    // -------------------------------------------------------------------------
    for (i = 0; i < n; i++)
    {
        SPEX_MPQ_SGN(&r, x->x.mpq[i]);
        if (r == 0) { continue;}

        for (int64_t p = 0; p < A->v[i]->nz; p++)
        {
            int64_t j = A->v[i]->i[p];
            SPEX_MPQ_SET_Z(temp, A->v[i]->x[p]);
            SPEX_MPQ_MUL(temp, temp, A->v[i]->scale);
            // b2[j] += x[i]*A(j,i)
            SPEX_MPQ_MUL(temp, temp, x->x.mpq[i]);
            SPEX_MPQ_ADD(b2->x.mpq[j], b2->x.mpq[j], temp);
        }
    }
    //--------------------------------------------------------------------------
    // Apply scales of A and b to b2 before comparing the b2 with scaled b'
    //--------------------------------------------------------------------------
    SPEX_MPQ_DIV(temp, b->scale, A->scale);

    // Apply scaling factor, but ONLY if it is not 1
    SPEX_MPQ_CMP_UI(&r, temp, 1, 1);
    if (r != 0)
    {
        for (i = 0; i < n; i++)
        {
            SPEX_MPQ_MUL(b2->x.mpq[i], b2->x.mpq[i], temp);
        }
    }

    // -------------------------------------------------------------------------
    // check if b2 == b
    // -------------------------------------------------------------------------
    *Is_correct = true;
    for (i = 0; i < n; i++)
    {
        // temp = b[i] (correct b)
        SPEX_MPQ_SET_Z(temp, b->x.mpz[i]);

        // set check false if b!=b2
        SPEX_MPQ_EQUAL(&r, temp, b2->x.mpq[i]);
        if (r == 0)
        {
            *Is_correct = false;
            break;
        }
    }

    //--------------------------------------------------------------------------
    // Print info
    //--------------------------------------------------------------------------

    if (*Is_correct)
    {
        printf ("Factorization is verified to be correct and exact.\n");
    }
    else
    {
        printf ("ERROR! Factorization is wrong. This is a bug; please "
                  "contact the authors of SPEX.\n");
    }

    SPEX_FREE_ALL;
    return info;
}
#endif
