//------------------------------------------------------------------------------
// SPEX_Utilities/SPEX_factorization_check: check the correctness of a given
// factorization.
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Chris Lourenco
// (US Naval Academy), Erick Moreno-Centeno, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// SPEX_factorization_check checks all the followings for a give factorization:
// 1. if the required components exist;
// 2. if sizes of different matrices match, i.e., L and U should be n*n, and
//    rhos should be n*1;
// 3. if L (and U if exists) is correct (using SPEX_matrix_check);
// 4. if L, U, and rhos have same pivot values, and when F is updatable, if L
//    (and U if exists) is of SPEX_DYNAMIC_CSC, and if the i-th pivot is the
//    first entry in the nonzero list of i-th vector of L (and U if exists);
// 5. if each permutation is reasonable, i.e., no duplicate, and in range of
//    [0,n), and if P_perm and Pinv_perm are mutually inverse vectors,
//    same applied to (Q_perm, Qinv_perm) if exists.

#define SPEX_FREE_ALL    \
{                        \
    SPEX_MPZ_CLEAR(tmpz);\
    SPEX_FREE(work);     \
}

#include "spex_util_internal.h"

// if pr == 2, turn off printing after 30 lines of output
#define SPEX_PR_LIMIT                       \
    lines++ ;                               \
    if (pr == 2 && lines > 30)              \
    {                                       \
        SPEX_PRINTF ("    ...\n") ;         \
        pr = 1 ;                            \
    }

SPEX_info SPEX_factorization_check
(
    SPEX_factorization F,       // The factorization to check
    const SPEX_options option
)
{

    //--------------------------------------------------------------------------
    // 1. basic check of the factorization object
    //--------------------------------------------------------------------------

    SPEX_info info = spex_factorization_basic_check (F) ;
    if (info != SPEX_OK) return (info) ;

    int pr = SPEX_OPTION_PRINT_LEVEL(option);

    //--------------------------------------------------------------------------
    // 2. check if sizes of different matrices match, i.e., L and U should be
    //    n*n, and rhos should be n*1;
    //--------------------------------------------------------------------------

    int64_t n = F->L->n;
    if (F->L->m != n || // L and U must be n-by-n
        (F->kind == SPEX_LU_FACTORIZATION && (F->U->m != n || F->U->n != n)) ||
        F->rhos->m != n || F->rhos->n != 1) // rhos must be n-by-1
    {
        return SPEX_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // 3. check if L, rhos (and U if exists) is correct
    //--------------------------------------------------------------------------

    int *work = NULL;
    mpz_t tmpz;
    SPEX_MPZ_SET_NULL(tmpz);

    SPEX_PR2("L:\n");
    SPEX_CHECK(SPEX_matrix_check (F->L, option)) ;

    if (F->kind == SPEX_LU_FACTORIZATION)
    {
        SPEX_PR2("U:\n");
        SPEX_CHECK(SPEX_matrix_check (F->U, option)) ;
    }

    SPEX_PR2("rhos:\n");
    SPEX_CHECK(SPEX_matrix_check (F->rhos, option)) ;

    //--------------------------------------------------------------------------
    // 4. check if each permutation is reasonable, i.e., no duplicate, and in
    //    range of [0,n), and if P_perm and Pinv_perm are mutually inverse
    //    vectors, same applied to (Q_perm, Qinv_perm) if exists.
    //--------------------------------------------------------------------------

    int64_t i, j, p;
    int r;
    int64_t lines = 0;     // # of lines printed so far

    // allocate workspace to check for duplicates
    work = (int *) SPEX_calloc (n, sizeof (int)) ;
    if (work == NULL)
    {
        // out of memory
        SPEX_PR1 ("out of memory\n") ;
        return (SPEX_OUT_OF_MEMORY) ;
    }
    for (j = 0; j < n; j++)
    {
        i = F->P_perm[j];
        if (i < 0 || i >= n)
        {
            // row indices out of range
            SPEX_PR1 ("index out of range: P(%ld) = %ld\n", j, i) ;
            SPEX_FREE_ALL ;
            return (SPEX_INCORRECT_INPUT) ;
        }
        else if (work [i] == 1)
        {
            // duplicate
            SPEX_PR1 ("duplicate index: P(%ld) = %ld\n", j, i) ;
            SPEX_FREE_ALL ;
            return (SPEX_INCORRECT_INPUT) ;
        }
        else if (F->Pinv_perm[i] != j)
        {
            // P_perm and Pinv_perm not mutually inverse
            SPEX_PR1 ("unmatched P and Pinv: Pinv(P(%ld)) = %ld\n", j,
                F->Pinv_perm[i]) ;
            SPEX_FREE_ALL ;
            return (SPEX_INCORRECT_INPUT) ;
        }
        work[i] = 1;
        SPEX_PR_LIMIT;
        SPEX_PR2("P[%ld] = %ld\n", j, i);
    }
    if (F->kind == SPEX_LU_FACTORIZATION)
    {
        lines = 0;
        for (j = 0; j < n; j++)
        {
            i = F->Q_perm[j];
            if (i < 0 || i >= n)
            {
                // row indices out of range
                SPEX_PR1 ("index out of range: Q(%ld) = %ld\n", j, i) ;
                SPEX_FREE_ALL ;
                return (SPEX_INCORRECT_INPUT) ;
            }
            else if (work [i] == 2)
            {
                // duplicate
                SPEX_PR1 ("duplicate index: Q(%ld) = %ld\n", j, i) ;
                SPEX_FREE_ALL ;
                return (SPEX_INCORRECT_INPUT) ;
            }
            else if (F->updatable && F->Qinv_perm[i] != j)
            {
                // Q_perm and Qinv_perm not mutually inverse
                SPEX_PR1 ("unmatched Q and Qinv: Qinv(Q(%ld)) = %ld\n", j,
                    F->Qinv_perm[i]) ;
                SPEX_FREE_ALL ;
                return (SPEX_INCORRECT_INPUT) ;
            }
            work[i] = 2;
            SPEX_PR_LIMIT;
            SPEX_PR2("Q[%ld] = %ld\n", j, i);
        }
    }

    //--------------------------------------------------------------------------
    // 5. check if L, U, and rhos have same pivot values, and when F is
    //    updatable, if L (and U if exists) is of SPEX_DYNAMIC_CSC, and if the
    //    i-th pivot is the first entry in the nonzero list of i-th vector of L
    //    (and U if exists);
    //--------------------------------------------------------------------------

    if (F->updatable)
    {
        SPEX_CHECK(SPEX_mpz_init(tmpz));
        for (i = 0; i < n; i++)
        {
            // first entry in the nnz pattern of i-th column of L must be pivot
            if (F->L->v[i]->i[0] != F->P_perm[i])
            {
                SPEX_PR1("first entry in L->v[%ld] is not pivot\n", i);
                SPEX_FREE_ALL;
                return SPEX_INCORRECT_INPUT;
            }
            else
            {
                // i-th pivot of L should be the same as rhos[i]
                SPEX_CHECK(SPEX_mpz_divexact(tmpz, F->L->v[i]->x[0],
                    SPEX_MPQ_DEN(F->L->v[i]->scale)));
                SPEX_CHECK(SPEX_mpz_mul(tmpz, tmpz,
                    SPEX_MPQ_NUM(F->L->v[i]->scale)));
                SPEX_CHECK(SPEX_mpz_cmp(&r, tmpz, F->rhos->x.mpz[i]));
                if (r != 0)
                {
                    SPEX_PR1("incorrect pivot: L(%ld, %ld) != rhos(%ld) \n",
                        i, i, i);
                    SPEX_FREE_ALL;
                    return SPEX_INCORRECT_INPUT;
                }
            }

            // for LU factorization, check for U as well
            if (F->kind == SPEX_LU_FACTORIZATION)
            {
                // first entry in the nnz pattern of i-th row of U must be pivot
                if (F->U->v[i]->i[0] != F->Q_perm[i])
                {
                    SPEX_PR1("first entry in U->v[%ld] is not pivot\n", i);
                    SPEX_FREE_ALL;
                    return SPEX_INCORRECT_INPUT;
                }
                else
                {
                    // i-th pivot of U should be the same as rhos[i]
                    SPEX_CHECK(SPEX_mpz_divexact(tmpz, F->U->v[i]->x[0],
                        SPEX_MPQ_DEN(F->U->v[i]->scale)));
                    SPEX_CHECK(SPEX_mpz_mul(tmpz, tmpz,
                        SPEX_MPQ_NUM(F->U->v[i]->scale)));
                    SPEX_CHECK(SPEX_mpz_cmp(&r, tmpz, F->rhos->x.mpz[i]));
                    if (r != 0)
                    {
                        SPEX_PR1("incorrect pivot: U(%ld, %ld) != rhos(%ld) \n",
                            i, i, i);
                        SPEX_FREE_ALL;
                        return SPEX_INCORRECT_INPUT;
                    }
                }
            }
        }
    }
    else // non-updatable factorization
    {
        for (i = 0; i < n; i++)
        {
            // get the index of the i-th pivot in L
            for (p = F->L->p[i]; p < F->L->p[i+1]; p++)
            {
                if (F->L->i[p] == i) break;
            }
            if (F->L->i[p] != i)
            {
                SPEX_PR1("L(%ld %ld) not found\n", i, i);
                SPEX_FREE_ALL;
                return SPEX_INCORRECT_INPUT;
            }

            // i-th pivot of L should be the same as rhos[i]
            SPEX_CHECK(SPEX_mpz_cmp(&r, F->L->x.mpz[p], F->rhos->x.mpz[i]));
            if (r != 0)
            {
                SPEX_PR1("incorrect pivot: L(%ld, %ld) != rhos(%ld) \n",
                    i, i, i);
                SPEX_FREE_ALL;
                return SPEX_INCORRECT_INPUT;
            }

            // for LU factorization, check for U as well
            if (F->kind == SPEX_LU_FACTORIZATION)
            {
                // get the index of the i-th pivot in U
                for (p = F->U->p[i]; p < F->U->p[i+1]; p++)
                {
                    if (F->U->i[p] == i) break;
                }
                if (F->U->i[p] != i)
                {
                    SPEX_PR1("U(%ld %ld) not found\n", i, i);
                    SPEX_FREE_ALL;
                    return SPEX_INCORRECT_INPUT;
                }

                // i-th pivot of U should be the same as rhos[i]
                SPEX_CHECK(SPEX_mpz_cmp(&r, F->U->x.mpz[p], F->rhos->x.mpz[i]));
                if (r != 0)
                {
                    SPEX_PR1("incorrect pivot: U(%ld, %ld) != rhos(%ld) \n",
                        i, i, i);
                    SPEX_FREE_ALL;
                    return SPEX_INCORRECT_INPUT;
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // All tests passed, free space and return SPEX_OK
    //--------------------------------------------------------------------------
    SPEX_FREE_ALL;
    return (SPEX_OK) ;
}

