//------------------------------------------------------------------------------
// SPEX_Update/spex_update_matrix_convert.c: convert between updatable and
// non-updatable matrix.
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Chris Lourenco (US Naval Academy),
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// spex_update_matrix_convert converts either F->L or F->U between CSC MPZ
// matrix and DYNAMIC_CSC MPZ matrix. This function is only called by
// spex_update_factorization_convert, and the updatable flag of the input
// factorization is always set to its complement, which indicates the expected
// result of the conversion. To be more specific, the result of the function is
// listed as follows:
// if F->updatable = true and convertL = true. This function converts F->L from
// CSC MPZ matrix to Dynamic_CSC MPZ matrix in an updatable format, which means
// that rows of L will be properly permuted with F->P_perm and the first entry
// of each column of L will be the corresponding digonal. 
// if F->updatable = true and convertL = false. This function converts F->U
// from CSC MPZ matrix to Dynamic_CSC MPZ matrix in an updatable format, which
// means that U will be firstly transposed to UT, and then the rows of UT will
// be properly permuted  with F->Q_perm, and the first entry of
// each column of UT will be the corresponding diagonal.
// if F->updatable = false and convertL = true. This function converts F->L
// from Dynamic_CSC MPZ matrix to CSC MPZ matrix in a non-updatable format,
// which means that the permutation of the rows of L will be reset with
// F->Pinv_perm.
// if F->updatable = false and convertL = false. This function converts F->U
// from Dynamic_CSC MPZ matrix to CSC MPZ matrix in a non-updatable format,
// which means that U will be firstly transposed to UT, and then the
// permutation of the rows of UT will be reset with F->Qinv_perm.
//
// If the function returns with any error, 

#define SPEX_FREE_ALL                \
{                                    \
    SPEX_FREE(rowcount);             \
    SPEX_FREE(Mp);                   \
    SPEX_FREE(Mi);                   \
    spex_delete_mpz_array(&Mx, nnz); \
    if (Mv != NULL)                  \
    {                                \
        for (i = 0; i < n; i++)      \
        {                            \
            SPEX_vector_free(&(Mv[i]), option);\
        }                            \
        SPEX_FREE(Mv);               \
    }                                \
}

#include "spex_update_internal.h"

SPEX_info spex_update_matrix_convert
(
    SPEX_factorization *F,       // converted CSC matrix
                                 // identity matrix if input as NULL
    const bool convertL,         // true if B->v[i] is the i-th col of B.
                                 // Otherwise, B->v[i] is the i-th row of B
    const SPEX_options *option
)
{
    SPEX_info info;
    SPEX_matrix *L = F->L, *U = F->U;
    SPEX_vector **Mv = NULL;
    mpz_t *Mx = NULL;
    int64_t i, j, p, m = L->m, n = L->n, mp, nnz = 0;
    int64_t *rowcount = NULL, *Mp = NULL, *Mi = NULL;
    bool found_diag = false;

    if (F->updatable)
    {
        //----------------------------------------------------------------------
        // allocate an array of n vetors to construct the DYNAMIC_CSC MPZ matrix
        //----------------------------------------------------------------------
        Mv = (SPEX_vector**) SPEX_calloc(n, sizeof(SPEX_vector*)); 
        if (!Mv) 
        { 
            SPEX_FREE_ALL ;
            return SPEX_OUT_OF_MEMORY; 
        } 
         
        for (int64_t i = 0; i < n; i++) 
        {
            SPEX_CHECK(SPEX_vector_allocate(&(Mv[i]), 0, option)); 
        } 
        if (convertL)
        {
            int64_t *Lp = L->p, *perm = F->P_perm;
            //------------------------------------------------------------------
            // convert L from CSC MPZ to DYNAMIC_CSC MPZ
            //------------------------------------------------------------------
            for (j = 0; j < n; j++)
            {
                // reallocate space for each column vector Mv[j]
                SPEX_CHECK(SPEX_vector_realloc(Mv[j], Lp[j+1]-Lp[j], option));
                mp = 0;
                found_diag = false;
                for (p = Lp[j]; p < Lp[j+1]; p++)
                {
                    i = L->i[p];
                    // permute entry indices
                    Mv[j]->i[mp] = perm[i];
                    // Mv[j]->x[mp] = L->x[p]
                    SPEX_CHECK(SPEX_mpz_swap(Mv[j]->x[mp], L->x.mpz[p]));

                    // make sure the first entry is the diagonal
                    if (!found_diag && i == j)
                    {
                        found_diag = true;
                        if (mp != 0)
                        {
                            SPEX_CHECK(SPEX_mpz_swap(Mv[j]->x[0],
                                                     Mv[j]->x[mp]));
                            Mv[j]->i[mp] = Mv[j]->i[0];
                            Mv[j]->i[0] = perm[i];
                        }
                    }

                    mp++;
                }
                Mv[j]->nz = mp;
            }

            // update components of L
            L->kind = SPEX_DYNAMIC_CSC;// type remains SPEX_MPZ
            SPEX_FREE(L->p);
            SPEX_FREE(L->i);
            spex_delete_mpz_array(&(L->x.mpz), L->nzmax);
            L->nzmax = 0;// reset nzmax
            L->v = Mv; // replace with Mv
            Mv = NULL;
        }
        else
        {
            int64_t *Up = U->p, *perm = F->Q_perm;
            //------------------------------------------------------------------
            // transpose U and convert it from CSC MPZ to DYNAMIC_CSC MPZ
            //------------------------------------------------------------------
            // find the # of nz in each row of U
            rowcount  = (int64_t*) SPEX_calloc(m, sizeof(int64_t));
            if (!rowcount)
            {
                SPEX_FREE_ALL;
                return SPEX_OUT_OF_MEMORY;
            }
            for (p = 0; p < Up[n]; p++)
            {
                rowcount[ U->i[p] ]++ ;
            }

            // reallocate space for each column vector Mv[i]
            for (i = 0; i < m; i++)
            {
                SPEX_CHECK(SPEX_vector_realloc(Mv[i], rowcount[i], option));
            }

            // construct Mv from U
            for (j = 0; j < n; j++)
            {
                found_diag = false;
                for (p = Up[j]; p < Up[j+1]; p++)
                {
                    i = U->i[p];
                    mp = Mv[i]->nz++;
                    Mv[i]->i[mp] = perm[j] ;
                    // Mv[i]->x[mp] = U->x[p]
                    SPEX_CHECK(SPEX_mpz_swap(Mv[i]->x[mp], U->x.mpz[p]));

                    // make sure the first entry is the diagonal
                    if (!found_diag && i == j)
                    {
                        found_diag = true;
                        if (mp != 0)
                        {
                            SPEX_CHECK(SPEX_mpz_swap(Mv[i]->x[0],
                                                     Mv[i]->x[mp]));
                            Mv[i]->i[mp] = Mv[i]->i[0];
                            Mv[i]->i[0] = perm[j];
                        }
                    }
                }
            }

            // update components of U
            U->kind = SPEX_DYNAMIC_CSC;// type remains SPEX_MPZ
            SPEX_FREE(U->p);
            SPEX_FREE(U->i);
            spex_delete_mpz_array(&(U->x.mpz), U->nzmax);
            U->nzmax = 0;// reset nzmax
            U->v = Mv; // replace with Mv
            Mv = NULL;
        }
    }
    else // convert to non-updatable format (CSC MPZ)
    {
        int sgn;
        if (convertL)
        {
            int64_t *perm = F->Pinv_perm;
            //------------------------------------------------------------------
            // convert L from DYNAMIC_CSC MPZ to CSC MPZ
            //------------------------------------------------------------------
            // compute number of nonzero in L
            for (j = 0; j < n; j++)
            {
                nnz += L->v[j]->nz;
            }

            // allocate space for Mi, Mp and Mx
            Mp = (int64_t *) SPEX_calloc (n+1, sizeof (int64_t)) ;
            Mi = (int64_t *) SPEX_calloc (nnz, sizeof (int64_t)) ;
            Mx = spex_create_mpz_array (nnz) ;
            if (!Mp || !Mi || !Mx)
            {
                SPEX_FREE_ALL;
                return SPEX_OUT_OF_MEMORY;
            }

            // construct Mi, Mp and Mx from L
            mp = 0;
            Mp[0] = 0;
            for (j = 0 ; j < n ; j++)
            {
                SPEX_CHECK(SPEX_mpq_cmp_ui(&sgn, L->v[j]->scale, 1, 1));
                for (p = 0 ; p < L->v[j]->nz ; p++)
                {
                    i = L->v[j]->i[p];
                    Mi[mp] = perm[i] ;
                    if (sgn != 0) // scale != 1
                    {
                        // apply scale to L->v[j]->x[p]
                        SPEX_CHECK(SPEX_mpz_divexact(L->v[j]->x[p],
                            L->v[j]->x[p], SPEX_MPQ_DEN(L->v[j]->scale))) ;
                        SPEX_CHECK(SPEX_mpz_mul(L->v[j]->x[p],
                            L->v[j]->x[p], SPEX_MPQ_NUM(L->v[j]->scale))) ;
                    }
                    SPEX_CHECK(SPEX_mpz_swap(Mx[mp], L->v[j]->x[p])) ;
                    mp++;
                }
                Mp[j+1] = mp;
            }

            // update components of L
            L->kind = SPEX_CSC;// type remains SPEX_MPZ
            L->p = Mp;     Mp = NULL;
            L->i = Mi;     Mi = NULL;
            L->x.mpz = Mx; Mx = NULL;
            L->nzmax = nnz;// update nzmax
            for (i = 0; i < n; i++)
            {
                SPEX_vector_free(&(L->v[i]), option);
            }
            SPEX_FREE(L->v);
        }
        else
        {
            int64_t *perm = F->Qinv_perm;
            //------------------------------------------------------------------
            // transpose U and convert it from DYNAMIC_CSC MPZ to CSC MPZ
            //------------------------------------------------------------------
            rowcount  = (int64_t*) SPEX_calloc(m, sizeof(int64_t));
            if (!rowcount)
            {
                SPEX_FREE_ALL;
                return SPEX_OUT_OF_MEMORY;
            }

            // find the # of nz in each row of U and total number of nonzero
            for (i = 0; i < n; i++)
            {
                nnz += U->v[i]->nz;
                for (p = 0 ; p < U->v[i]->nz; p++)
                {
                    j = U->v[i]->i[p];
                    rowcount [j]++ ;
                }
            }

            // allocate space for Mi, Mp and Mx
            Mp = (int64_t *) SPEX_calloc (n+1, sizeof (int64_t)) ;
            Mi = (int64_t *) SPEX_calloc (nnz, sizeof (int64_t)) ;
            Mx = spex_create_mpz_array (nnz) ;
            if (!Mp || !Mi || !Mx)
            {
                SPEX_FREE_ALL;
                return SPEX_OUT_OF_MEMORY;
            }

            // compute cumulative sum of rowcount to get the col pointer
            SPEX_CHECK(SPEX_cumsum(Mp, rowcount, m, option));

            // construct Mi, Mp and Mx from U
            for (i = 0 ; i < n ; i++)
            {
                SPEX_CHECK(SPEX_mpq_cmp_ui(&sgn, U->v[i]->scale, 1, 1));
                for (p = 0 ; p < U->v[i]->nz ; p++)
                {
                    j = U->v[i]->i[p];
                    mp = rowcount[j];
                    Mi[ mp ] = perm[i] ;
                    if (sgn != 0) // scale != 1
                    {
                        // apply scale to U->v[i]->x[p]
                        SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[p],
                            U->v[i]->x[p], SPEX_MPQ_DEN(U->v[i]->scale))) ;
                        SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[p],
                            U->v[i]->x[p], SPEX_MPQ_NUM(U->v[i]->scale))) ;
                    }
                    SPEX_CHECK(SPEX_mpz_swap(Mx[mp], U->v[i]->x[p])) ;
                    rowcount[j] ++;
                }
            }

            // update components of U
            U->kind = SPEX_CSC;// type remains SPEX_MPZ
            U->p = Mp;     Mp = NULL;
            U->i = Mi;     Mi = NULL;
            U->x.mpz = Mx; Mx = NULL;
            U->nzmax = nnz;// update nzmax
            for (i = 0; i < n; i++)
            {
                SPEX_vector_free(&(U->v[i]), option);
            }
            SPEX_FREE(U->v);
        }
    }

    SPEX_FREE_ALL;
    return SPEX_OK;
}
