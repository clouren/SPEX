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

#define SPEX_FREE_ALL                \
    SPEX_FREE(rowcount);             \
    SPEX_matrix_free(&M, option);

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
    SPEX_matrix *M = NULL, *L = F->L, *U = F->U;
    int64_t i, j, p, m = L->m, n = L->n, Mp;
    int64_t *rowcount = NULL;
    bool found_diag = false;

    if (F->updatable)
    {
        //----------------------------------------------------------------------
        // allocate M as DYNAMIC_CSC MPZ matrix, which contains n empty vectors
        //----------------------------------------------------------------------
        SPEX_CHECK(SPEX_matrix_allocate(&M, SPEX_DYNAMIC_CSC, SPEX_MPZ, m, n,
            0, false, true, option));
        if (convertL)
        {
            int64_t *Lp = L->p, *perm = F->P_perm;
            //------------------------------------------------------------------
            // convert L from CSC MPZ to DYNAMIC_CSC MPZ
            //------------------------------------------------------------------
            for (j = 0; j < n; j++)
            {
                // reallocate space for each column of A
                SPEX_CHECK(SPEX_vector_realloc(M->v[j], Lp[j+1]-Lp[j], option));
                Mp = 0;
                found_diag = false;
                for (p = Lp[j]; p < Lp[j+1]; p++)
                {
                    i = L->i[p];
                    // permute entry indices
                    M->v[j]->i[Mp] = perm[i];
                    // M->v[j]->x[Mp] = L->x[p]
                    SPEX_CHECK(SPEX_mpz_swap(M->v[j]->x[Mp], L->x.mpz[p]));

                    // make sure the first entry is the diagonal
                    if (!found_diag && i == j)
                    {
                        found_diag = true;
                        if (Mp != 0)
                        {
                            SPEX_CHECK(SPEX_mpz_swap(M->v[j]->x[0],
                                                     M->v[j]->x[Mp]));
                            M->v[j]->i[Mp] = M->v[j]->i[0];
                            M->v[j]->i[0] = perm[i];
                        }
                    }

                    Mp++;
                }
                M->v[j]->nz = Mp;
            }

            // replace original L matrix with M
            SPEX_CHECK(SPEX_matrix_free(&(F->L), option));
            F->L = M;
            M = NULL;
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

            // reallocate space for each column of M
            for (i = 0; i < m; i++)
            {
                SPEX_CHECK(SPEX_vector_realloc(M->v[i], rowcount[i], option));
            }

            // construct M from U
            for (j = 0; j < n; j++)
            {
                found_diag = false;
                for (p = Up[j]; p < Up[j+1]; p++)
                {
                    i = U->i[p];
                    Mp = M->v[i]->nz++;
                    M->v[i]->i[Mp] = perm[j] ;
                    // M->v[i]->x[Mp] = U->x[p]
                    SPEX_CHECK(SPEX_mpz_swap(M->v[i]->x[Mp], U->x.mpz[p]));

                    // make sure the first entry is the diagonal
                    if (!found_diag && i == j)
                    {
                        found_diag = true;
                        if (Mp != 0)
                        {
                            SPEX_CHECK(SPEX_mpz_swap(M->v[i]->x[0],
                                                     M->v[i]->x[Mp]));
                            M->v[i]->i[Mp] = M->v[i]->i[0];
                            M->v[i]->i[0] = perm[j];
                        }
                    }
                }
            }
            // replace original U matrix with M
            SPEX_CHECK(SPEX_matrix_free(&(F->U), option));
            F->U = M;
            M = NULL;
        }
    }
    else // convert to non-updatable format (CSC MPZ)
    {
        int64_t nnz = 0;
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

            // allocate space for M
            SPEX_CHECK(SPEX_matrix_allocate(&M, SPEX_CSC, SPEX_MPZ, m, n, nnz,
                false, true, option));

            // initialize for M and construct M from L
            nnz = 0;
            M->p[0] = 0;
            for (j = 0 ; j < n ; j++)
            {
                SPEX_CHECK(SPEX_mpq_cmp_ui(&sgn, L->v[j]->scale, 1, 1));
                for (p = 0 ; p < L->v[j]->nz ; p++)
                {
                    i = L->v[j]->i[p];
                    M->i[nnz] = perm[i] ;
                    if (sgn != 0) // scale != 1
                    {
                        // apply scale to L->v[j]->x[p]
                        SPEX_CHECK(SPEX_mpz_divexact(L->v[j]->x[p],
                            L->v[j]->x[p], SPEX_MPQ_DEN(L->v[j]->scale))) ;
                        SPEX_CHECK(SPEX_mpz_mul(L->v[j]->x[p],
                            L->v[j]->x[p], SPEX_MPQ_NUM(L->v[j]->scale))) ;
                    }
                    SPEX_CHECK(SPEX_mpz_swap(M->x.mpz[nnz], L->v[j]->x[p])) ;
                    nnz++;
                }
                M->p[j+1] = nnz;
            }

            // replace original L matrix with M
            SPEX_CHECK(SPEX_matrix_free(&(F->L), option));
            F->L = M;
            M = NULL;
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

            // allocate space for M
            SPEX_CHECK(SPEX_matrix_allocate(&M, SPEX_CSC, SPEX_MPZ, m, n, nnz,
                false, true, option));

            // compute cumulative sum of rowcount to get the col pointer
            SPEX_CHECK(SPEX_cumsum(M->p, rowcount, m));

            // construct M from U
            for (i = 0 ; i < n ; i++)
            {
                SPEX_CHECK(SPEX_mpq_cmp_ui(&sgn, U->v[i]->scale, 1, 1));
                for (p = 0 ; p < U->v[i]->nz ; p++)
                {
                    j = U->v[i]->i[p];
                    Mp = rowcount[j];
                    M->i[ Mp ] = perm[i] ;
                    if (sgn != 0) // scale != 1
                    {
                        // apply scale to U->v[i]->x[p]
                        SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[p],
                            U->v[i]->x[p], SPEX_MPQ_DEN(U->v[i]->scale))) ;
                        SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[p],
                            U->v[i]->x[p], SPEX_MPQ_NUM(U->v[i]->scale))) ;
                    }
                    SPEX_CHECK(SPEX_mpz_swap(M->x.mpz[Mp], U->v[i]->x[p])) ;
                    rowcount[j] ++;
                }
            }

            // replace original U matrix with M
            SPEX_CHECK(SPEX_matrix_free(&(F->U), option));
            F->U = M;
            M = NULL;
        }
    }

    SPEX_FREE_ALL;
    return SPEX_OK;
}
