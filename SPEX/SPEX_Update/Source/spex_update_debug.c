//------------------------------------------------------------------------------
// SPEX_Update/spex_update_debug: for debugging purpose
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2022, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Purpose: This function performs direct LU factorization for A, and use
// this direct factorization as base line to verify if input factorization


#define SPEX_FREE_ALL                          \
    SPEX_FREE(opt1);                           \
    SPEX_symbolic_analysis_free(&Stmp, option);\
    SPEX_matrix_free(&Atmp, option);           \
    SPEX_matrix_free(&UT, option);             \
    SPEX_MPQ_CLEAR(tmpq);                      \
    SPEX_MPZ_CLEAR(tmpz);                      \
    SPEX_factorization_free(&Ftmp, option);

#include "spex_update_internal.h"

#ifdef SPEX_DEBUG
// This function is only compiled in debug mode.

#define SL(k) (F->L->v[(k)]->scale)
#define SU(k) (F->U->v[(k)]->scale)

SPEX_info spex_update_debug
(
    bool *Is_correct,       // if factorization is correct
    SPEX_factorization F,   // LU factorization of A
    int64_t k,              // current iteration
    spex_scattered_vector Lk_dense_col, // scattered column k of L
    spex_scattered_vector Uk_dense_row, // scattered column k of U
    bool finish_update,     // if the update process has finished
    SPEX_matrix vk,         // inserted column
    SPEX_matrix A,          // Input matrix Dynamic_CSC MPZ
    const SPEX_options option   // Command parameters
)
{

    SPEX_info info = SPEX_OK;
    *Is_correct = true;
    int64_t n = F->L->n, *P = F->P_perm, *Q = F->Q_perm;
    int r;
    SPEX_symbolic_analysis Stmp = NULL;
    SPEX_matrix Atmp = NULL;
    SPEX_matrix UT = NULL;
    SPEX_factorization Ftmp = NULL;
    SPEX_options opt1 = NULL;
    mpq_t tmpq; SPEX_MPQ_SET_NULL(tmpq);
    mpz_t tmpz; SPEX_MPZ_SET_NULL(tmpz);

    Stmp = (SPEX_symbolic_analysis) SPEX_calloc(1,
        sizeof(SPEX_symbolic_analysis_struct));
    if (Stmp == NULL) {return SPEX_OUT_OF_MEMORY;}
    Stmp->kind = SPEX_LU_FACTORIZATION;

    //--------------------------------------------------------------------------
    // create symbolic analysis
    //--------------------------------------------------------------------------
    // Allocate memory for column permutation
    Stmp->Q_perm = (int64_t*) SPEX_malloc((n+1) * sizeof(int64_t));
    if (Stmp->Q_perm == NULL)
    {
        SPEX_FREE(Stmp);
        SPEX_CHECK(SPEX_OUT_OF_MEMORY);
    }
    memcpy(Stmp->Q_perm, F->Q_perm, n*sizeof(int64_t));
    SPEX_CHECK(SPEX_matrix_nnz(&(Stmp->lnz), F->L, option));
    SPEX_CHECK(SPEX_matrix_nnz(&(Stmp->unz), F->U, option));

    //--------------------------------------------------------------------------
    // copy and permute A
    //--------------------------------------------------------------------------
    if(finish_update)
    {
        SPEX_CHECK(SPEX_update_matrix_colrep(A, vk, k, option));
        SPEX_CHECK(spex_update_verify(Is_correct, F, A, option));
        if (*Is_correct)
        {
            SPEX_CHECK(SPEX_update_matrix_colrep(A, vk, k, option));
            SPEX_FREE_ALL;
            return SPEX_OK;
        }
        SPEX_CHECK(SPEX_matrix_copy(&Atmp, SPEX_CSC, SPEX_MPZ, A, option));
        SPEX_CHECK(SPEX_update_matrix_colrep(A, vk, k, option));
    }
    else
    {
        SPEX_CHECK(SPEX_matrix_copy(&Atmp, SPEX_CSC, SPEX_MPZ, A, option));
    }
    for (int64_t jcol = 0; jcol < n; jcol++)
    {
        for (int64_t p = Atmp->p[jcol]; p < Atmp->p[jcol+1]; p++)
        {
            int64_t irow = Atmp->i[p];
            Atmp->i[p] = F->Q_perm[F->Pinv_perm[irow]];
        }
    }

    //--------------------------------------------------------------------------
    // perform LU factorization for permuted A with no ordering
    //--------------------------------------------------------------------------
    SPEX_CHECK(SPEX_create_default_options(&opt1));
    opt1->pivot = SPEX_DIAGONAL;
    SPEX_CHECK(SPEX_lu_factorize(&Ftmp, Atmp, Stmp, opt1));

    SPEX_MPQ_INIT(tmpq);
    for (int64_t jcol = 0; jcol < n; jcol++)
    {
        int64_t p, p1;
        bool fail = false;
        //----------------------------------------------------------------------
        // check for rhos
        //----------------------------------------------------------------------
        SPEX_MPZ_CMP(&r, Ftmp->rhos->x.mpz[jcol],
                                F->rhos->x.mpz[jcol]);
        if (r != 0)
        {
            mpq_set_num(tmpq, Ftmp->rhos->x.mpz[jcol]);
            mpq_set_den(tmpq, F->rhos->x.mpz[jcol]);
            mpq_canonicalize(tmpq);
            gmp_printf("correct/updated rhos ratio=%Qd\n",tmpq);
            printf("rhos(%ld) incorrect, ratio = %f\n", jcol,
                mpq_get_d(tmpq));
            *Is_correct = false;
            SPEX_FREE_ALL;
            return SPEX_OK;
        }

        //----------------------------------------------------------------------
        // check if pivots of L match rhos
        //----------------------------------------------------------------------
        if (jcol == k && !finish_update)
        {
            SPEX_MPZ_DIVEXACT(tmpz, Lk_dense_col->x[P[jcol]],
                                         SPEX_MPQ_DEN(SL(jcol)));
        }
        else
        {
            SPEX_MPZ_DIVEXACT(tmpz, F->L->v[jcol]->x[0],
                                         SPEX_MPQ_DEN(SL(jcol)));
        }
        SPEX_MPZ_MUL     (tmpz, tmpz,
                                     SPEX_MPQ_NUM(SL(jcol)));
        SPEX_MPZ_CMP(&r, Ftmp->rhos->x.mpz[jcol], tmpz);
        if (r != 0)
        {
            printf("pivot of %ld-th column of L does not match rhos\n", jcol);
            *Is_correct = false;
            SPEX_FREE_ALL;
            return SPEX_OK;
        }
        //----------------------------------------------------------------------
        // check if pivots of U match rhos
        //----------------------------------------------------------------------
        if (jcol == k && !finish_update)
        {
            SPEX_MPZ_DIVEXACT(tmpz, Uk_dense_row->x[Q[jcol]],
                                         SPEX_MPQ_DEN(SU(jcol)));
        }
        else
        {
            SPEX_MPZ_DIVEXACT(tmpz, F->U->v[jcol]->x[0],
                                         SPEX_MPQ_DEN(SU(jcol)));
        }
        SPEX_MPZ_MUL     (tmpz, tmpz,
                                     SPEX_MPQ_NUM(SU(jcol)));
        SPEX_MPZ_CMP(&r, Ftmp->rhos->x.mpz[jcol], tmpz);
        if (r != 0)
        {
            printf("pivot of %ld-th row of U does not match rhos\n", jcol);
            *Is_correct = false;
            SPEX_FREE_ALL;
            return SPEX_OK;
        }

        //----------------------------------------------------------------------
        // check for values of all entries in L
        //----------------------------------------------------------------------

        for (p = Ftmp->L->p[jcol]; p < Ftmp->L->p[jcol+1]; p++)
        {
            int64_t irow = Ftmp->L->i[p];
            int64_t found_p = -1;
            if (jcol == k && !finish_update)
            {
                for (p1 = 0; p1 < Lk_dense_col->nz; p1++)
                {
                    if (P[irow] == Lk_dense_col->i[p1])
                    {
                        found_p = p1;
                        break;
                    }
                }
            }
            else
            {
                for (p1 = 0; p1 < F->L->v[jcol]->nz; p1++)
                {
                    if (P[irow] == F->L->v[jcol]->i[p1])
                    {
                        found_p = p1;
                        break;
                    }
                }
            }
            if (found_p == -1)
            {
                SPEX_MPZ_SGN(&r, Ftmp->L->x.mpz[p]);
                if (r != 0)
                {
                    printf("L(%ld,%ld) should not be 0\n", irow, jcol);
                    fail = true;
                }
            }
            else
            {
                if (jcol == k && !finish_update)
                {
                    SPEX_MPZ_DIVEXACT(tmpz,
                                                 Lk_dense_col->x[P[irow]],
                                                 SPEX_MPQ_DEN(SL(jcol)));
                }
                else
                {
                    SPEX_MPZ_DIVEXACT(tmpz,
                                                 F->L->v[jcol]->x[found_p],
                                                 SPEX_MPQ_DEN(SL(jcol)));
                }
                SPEX_MPZ_MUL     (tmpz, tmpz,
                                             SPEX_MPQ_NUM(SL(jcol)));
                SPEX_MPZ_CMP(&r, Ftmp->L->x.mpz[p], tmpz);
                if (r != 0)
                {
                    mpq_set_num(tmpq, Ftmp->L->x.mpz[p]);
                    mpq_set_den(tmpq, tmpz);
                    mpq_canonicalize(tmpq);
                    printf("L(%ld,%ld) incorrect, ratio = %f\n", irow,
                        jcol, mpq_get_d(tmpq));
                    fail = true;
                }
            }
        }
        if (fail)
        {
            printf("expected %ld-th column of L:\n",jcol);
            for (p = Ftmp->L->p[jcol]; p < Ftmp->L->p[jcol+1]; p++)
            {
                if (mpz_sgn(Ftmp->L->x.mpz[p]) == 0) continue;
                int64_t irow = Ftmp->L->i[p];
                printf("%ld ",irow);
            }
            printf("\n--------------------\n");
            printf("updated %ld-th column of L:\n",jcol);
            if (jcol == k && !finish_update)
            {
                for (p = 0; p < Lk_dense_col->nz; p++)
                {
                    int64_t irow = Lk_dense_col->i[p];
                    if (mpz_sgn(Lk_dense_col->x[irow]) == 0) continue;
                    printf("%ld ",F->Pinv_perm[irow]);
                }
            }
            else
            {
                for (p = 0; p < F->L->v[jcol]->nz; p++)
                {
                    if (mpz_sgn(F->L->v[jcol]->x[p]) == 0) continue;
                    int64_t irow = F->L->v[jcol]->i[p];
                    printf("%ld ",F->Pinv_perm[irow]);
                }
            }
            printf("\n\n");
            *Is_correct = false;
            SPEX_FREE_ALL;
            return SPEX_OK;
        }
    }

    //----------------------------------------------------------------------
    // check for values of all entries in U
    //----------------------------------------------------------------------
    SPEX_CHECK(SPEX_transpose(&UT, Ftmp->U, option));

    for (int64_t i = 0; i < n; i++)
    {
        int64_t p, p1;
        bool fail = false;

        for (p = UT->p[i]; p < UT->p[i+1]; p++)
        {
            int64_t j = UT->i[p];
            int64_t found_p = -1;
            if (i == k && !finish_update)
            {
                for (p1 = 0; p1 < Uk_dense_row->nz; p1++)
                {
                    if (Q[j] == Uk_dense_row->i[p1])
                    {
                        found_p = p1;
                        break;
                    }
                }
            }
            else
            {
                for (p1 = 0; p1 < F->U->v[i]->nz; p1++)
                {
                    if (Q[j] == F->U->v[i]->i[p1])
                    {
                        found_p = p1;
                        break;
                    }
                }
            }
            if (found_p == -1)
            {
                SPEX_MPZ_SGN(&r, UT->x.mpz[p]);
                if (r != 0)
                {
                    // entries in the column k of U have been deleted
                    if (finish_update || j != k)
                    {
                        printf("U(%ld,%ld) should not be 0\n", i, j);
                        fail = true;
                    }
                }
            }
            else
            {
                if (i == k && !finish_update)
                {
                    SPEX_MPZ_DIVEXACT(tmpz, Uk_dense_row->x[Q[j]],
                                                 SPEX_MPQ_DEN(SU(i)));
                }
                else
                {
                    SPEX_MPZ_DIVEXACT(tmpz, F->U->v[i]->x[found_p],
                                                 SPEX_MPQ_DEN(SU(i)));
                }
                SPEX_MPZ_MUL     (tmpz, tmpz,
                                             SPEX_MPQ_NUM(SU(i)));
                SPEX_MPZ_CMP(&r, UT->x.mpz[p], tmpz);
                if (r != 0)
                {
                    mpq_set_num(tmpq, UT->x.mpz[p]);
                    mpq_set_den(tmpq, tmpz);
                    mpq_canonicalize(tmpq);
                    printf("U(%ld,%ld) incorrect, ratio = %f\n", i, j,
                        mpq_get_d(tmpq));
                    fail = true;
                }
            }
        }
        if (fail)
        {
            printf("expected %ld-th row of U:\n",i);
            for (p = UT->p[i]; p < UT->p[i+1]; p++)
            {
                if (mpz_sgn(UT->x.mpz[p]) == 0) continue;
                int64_t j = UT->i[p];
                printf("%ld ",j);
            }
            printf("\n--------------------\n");
            printf("updated %ld-th row of U:\n",i);
            if (i == k && !finish_update)
            {
                for (p = 0; p < Uk_dense_row->nz; p++)
                {
                    int64_t j = Uk_dense_row->i[p];
                    if (mpz_sgn(Uk_dense_row->x[j]) == 0) continue;
                    printf("%ld ",F->Qinv_perm[j]);
                }
            }
            else
            {
                for (p = 0; p < F->U->v[i]->nz; p++)
                {
                    if (mpz_sgn(F->U->v[i]->x[p]) == 0) continue;
                    int64_t j = F->U->v[i]->i[p];
                    printf("%ld ",F->Qinv_perm[j]);
                }
            }
            printf("\n\n");
            *Is_correct = false;
            SPEX_FREE_ALL;
            return SPEX_OK;
        }
    }

    SPEX_FREE_ALL;
    return info;
}

#endif
