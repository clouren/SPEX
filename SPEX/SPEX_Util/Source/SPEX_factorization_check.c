//------------------------------------------------------------------------------
// SPEX_Util/SPEX_factorization_check.c: check the correctness of a given
// factorization.
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Chris Lourenco
// (US Naval Academy), Erick Moreno-Centeno, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// SPEX_factorization_check checks all the followings for a give factorization:
// 1. if the required components exist;
//    (TODO use spex_factorization_basic_check?)
// 2. if sizes of different matrices match, i.e., L and U should be n*n, and
//    rhos should be n*1;
// 3. if L (and U if exists) is correct (using SPEX_matrix_check);
// 4. if L, U, and rhos have same pivot values, and when F is updatable, if L
//    (and U if exists) is of SPEX_DYNAMIC_CSC, and if the i-th pivot is the
//    first entry in the nonzero list of i-th vector of L (and U if exists);
// 5. if each permutation is reasonable, i.e., no duplicate, and in range of
//    [0,n), and if P_perm and Pinv_perm are mutually inverse vectors,
//    same applied to (Q_perm, Qinv_perm) if exists.

#define SPEX_FREE_ALL   \
{                       \
    SPEX_FREE(work);    \
}

#include "spex_util_internal.h"

SPEX_info SPEX_factorization_check
(
    SPEX_factorization *F, // The factorization to check
    const SPEX_options* option
)
{

    if (!spex_initialized()) {return SPEX_PANIC;}

    //--------------------------------------------------------------------------
    // 1. check if the required components exist;
    //--------------------------------------------------------------------------

    SPEX_info info;
    int pr = SPEX_OPTION_PRINT_LEVEL(option);

    // TODO: maybe put this if test in a helper function
    // spex_factorization_basic_check (F)
    if (!F ||
        /* F can only be LU or Cholesky*/
        !(F->kind == SPEX_LU_FACTORIZATION ||
          F->kind == SPEX_CHOLESKY_FACTORIZATION) ||
        /* rhos, P_perm and Pinv_perm should exist for any kind*/ 
        !(F->rhos) || !(F->P_perm) || !(F->Pinv_perm) ||
        /* L should exist and have MPZ entries*/
        !(F->L) || F->L->type != SPEX_MPZ ||
        /* L should be non-shallow csc for non-updatable or
                           dynamic csc for updatable*/
        !(
          (!(F->updatable) && F->L->kind == SPEX_CSC &&
                              !(F->L->p_shallow) && F->L->p &&
                              !(F->L->i_shallow) && F->L->i &&
                              !(F->L->x_shallow) && F->L->x.mpz     ) ||
          (  F->updatable  && F->L->kind == SPEX_DYNAMIC_CSC && F->L->v)
         ) ||
        /* if F is LU*/
        (F->kind == SPEX_LU_FACTORIZATION &&
        /* Q_perm should always exist and Qinv_perm should exist for updatable*/
         (!(F->Q_perm) || (F->updatable && !(F->Qinv_perm)) ||
        /*U should exist and have MPZ entries*/
          !(F->U) || F->U->type != SPEX_MPZ ||
        /* U should be non-shallow csc for non-updatable or
                           dynamic csc for updatable*/
          !(
            (!(F->updatable) && F->U->kind == SPEX_CSC &&
                                !(F->U->p_shallow) && F->U->p &&
                                !(F->U->i_shallow) && F->U->i &&
                                !(F->U->x_shallow) && F->U->x.mpz     ) ||
            (  F->updatable  && F->U->kind == SPEX_DYNAMIC_CSC && F->U->v)
           )                 )))
    {
        return SPEX_INCORRECT_INPUT;
    }


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
    SPEX_CHECK(SPEX_matrix_check (F->L, option)) ;
    SPEX_CHECK(SPEX_matrix_check (F->rhos, option)) ;

    if (F->kind == SPEX_LU_FACTORIZATION)
    {
        SPEX_CHECK(SPEX_matrix_check (F->U, option)) ;
    }

    //--------------------------------------------------------------------------
    // 4. check if L, U, and rhos have same pivot values, and when F is
    //    updatable, if L (and U if exists) is of SPEX_DYNAMIC_CSC, and if the
    //    i-th pivot is the first entry in the nonzero list of i-th vector of L
    //    (and U if exists);
    //--------------------------------------------------------------------------

    int64_t i, j, p;
    int r;
    char *buff = NULL;
    SPEX_info status = 0;
    if (F->updatable)
    {
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
                SPEX_CHECK(SPEX_mpz_cmp(&r, F->L->v[i]->x[0],
                                            F->rhos->x.mpz[i]));
                if (r != 0)
                {
                    if (pr >= 1)
                    {
                        // FIXME use gmp_printf
                        status = SPEX_mpfr_asprintf(&buff,
                            "L(%ld, %ld) = %Zd != rhos(%ld) = %Zd \n",
                            i, i, F->L->v[i]->x[0], i, F->rhos->x.mpz [i]);
                        if (status < 0)
                        {
                            SPEX_FREE_ALL ;
                            SPEX_PRINTF (" error: %d\n", status) ;
                            return (status) ;
                        }
                        else
                        {
                            SPEX_PR1("%s", buff);
                            SPEX_mpfr_free_str (buff);
                        }
                    }

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
                    SPEX_CHECK(SPEX_mpz_cmp(&r, F->U->v[i]->x[0],
                                                F->rhos->x.mpz[i]));
                    if (r != 0)
                    {
                        if (pr >= 1)
                        {
                            // FIXME use gmp_printf
                            status = SPEX_mpfr_asprintf(&buff,
                                "U(%ld, %ld) = %Zd != rhos(%ld) = %Zd \n",
                                i, i, F->U->v[i]->x[0], i, F->rhos->x.mpz [i]);
                            if (status < 0)
                            {
                                SPEX_FREE_ALL ;
                                SPEX_PRINTF (" error: %d\n", status) ;
                                return (status) ;
                            }
                            else
                            {
                                SPEX_PR1("%s", buff);
                                SPEX_mpfr_free_str (buff);
                            }
                        }
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
                if (pr >= 1)
                {
                    // FIXME use gmp_printf
                    status = SPEX_mpfr_asprintf(&buff,
                        "L(%ld, %ld) = %Zd != rhos(%ld) = %Zd \n",
                        i, i, F->L->x.mpz[p], i, F->rhos->x.mpz [i]);
                    if (status < 0)
                    {
                        SPEX_FREE_ALL ;
                        SPEX_PRINTF (" error: %d\n", status) ;
                        return (status) ;
                    }
                    else
                    {
                        SPEX_PR1("%s", buff);
                        SPEX_mpfr_free_str (buff);
                    }
                }
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
                    if (pr >= 1)
                    {
                        // FIXME use gmp_printf
                        status = SPEX_mpfr_asprintf(&buff,
                            "U(%ld, %ld) = %Zd != rhos(%ld) = %Zd \n",
                            i, i, F->U->x.mpz[p], i, F->rhos->x.mpz [i]);
                        if (status < 0)
                        {
                            SPEX_FREE_ALL ;
                            SPEX_PRINTF (" error: %d\n", status) ;
                            return (status) ;
                        }
                        else
                        {
                            SPEX_PR1("%s", buff);
                            SPEX_mpfr_free_str (buff);
                        }
                    }
                    SPEX_FREE_ALL;
                    return SPEX_INCORRECT_INPUT;
                }
            }
        }
    }

    //--------------------------------------------------------------------------
    // 5. check if each permutation is reasonable, i.e., no duplicate, and in
    //    range of [0,n), and if P_perm and Pinv_perm are mutually inverse
    //    vectors, same applied to (Q_perm, Qinv_perm) if exists.
    //--------------------------------------------------------------------------
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
    }
    if (F->kind == SPEX_LU_FACTORIZATION)
    {
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
        }
    }

    //--------------------------------------------------------------------------
    // All tests passed, free space and return SPEX_OK
    //--------------------------------------------------------------------------
    SPEX_FREE_ALL;
    return (SPEX_OK) ;
}

