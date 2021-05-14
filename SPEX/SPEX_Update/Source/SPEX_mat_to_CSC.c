//------------------------------------------------------------------------------
// SPEX_Update/SPEX_mat_to_CSC.c: convert SPEX_mat matrix to a CSC matrix
// stored as SPEX_matrix
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// SPEX_mat_to_CSC create a SPEX_matrix A in Compressed Sparse Column (CSC)
// format with entries in mpz_t type from the give matrix B, which is used in
// the SPEX LU update process. The entries of B are also in mpz_t type but the
// entries is stored either in column form or in row form, which is specified
// by B_Is_ColWise. Since the SPEX LU update process requires A = LDU, while
// the SPEX LU factorization gives PAQ = LDU, additional permutation will be
// applied.
//
// When converting A, we should use SPEX_mat_to_CSC(A, A1, NULL, true, option)
// and the resulting A will be in the same column and row order as A1.
//
// When converting L, we should use SPEX_mat_to_CSC(L, L1, P_inv, true, option)
// and the resulting L will be strictly a lower triangular matrix.
//
// When converting U, we should use SPEX_mat_to_CSC(U, U1, Q_inv, false, option)
// and the resulting U will be a strictly upper triangular matrix.

#define SPEX_FREE_WORK               \
    SPEX_FREE(rowcount);

#define SPEX_FREE_ALL                \
    SPEX_FREE_WORK;                  \
    SPEX_matrix_free(&A, option);

#include "spex_update_internal.h"

SPEX_info SPEX_mat_to_CSC
(
    SPEX_matrix **A_handle,       // converted CSC matrix
    const SPEX_mat *B,            // original matrix
    const int64_t *perm_inv,      // inverse of permutation, consider as
                                  // identity matrix if input as NULL
    const bool B_Is_ColWise,      // true if B->v[i] is the i-th col of B.
                                  // Otherwise, B->v[i] is the i-th row of B
    const SPEX_options *option
)
{
    if (A_handle == NULL || B == NULL)
    {
        return SPEX_INCORRECT_INPUT;
    }
    SPEX_info info;
    int sgn;
    int64_t *rowcount = NULL;
    (*A_handle) = NULL;
    SPEX_matrix *A = NULL;

    int64_t i, j, p, nnz = 0;

    if (B_Is_ColWise)
    {
        // if B is stored in a column-wise format, find # of nnz in B first
        for (j = 0; j < B->n; j++)
        {
            nnz += B->v[j]->nz;
        }

        // allocate space for A
        SPEX_CHECK(SPEX_matrix_allocate(&A, SPEX_CSC, SPEX_MPZ, B->n, B->m, nnz,
            false, true, option));

        // initialize for A and construct A from B
        nnz = 0;
        A->p[0] = 0;
        for (j = 0 ; j < B->n ; j++)
        {
            SPEX_CHECK(SPEX_mpq_cmp_ui(&sgn, B->v[j]->scale, 1, 1));
            for (p = 0 ; p < B->v[j]->nz ; p++)
            {
                i = B->v[j]->i[p];
                A->i[nnz] = (perm_inv == NULL) ? i : perm_inv[i] ;
                // A->x[nnz] = B->v[j]->x[p]*scale
                if (sgn == 0) // scale == 1
                {
                    SPEX_CHECK(SPEX_mpz_set(SPEX_1D(A, nnz, mpz),
                        B->v[j]->x[p])) ;
                }
                else
                {
                    SPEX_CHECK(SPEX_mpz_divexact(SPEX_1D(A, nnz, mpz),
                        B->v[j]->x[p], SPEX_MPQ_DEN(B->v[j]->scale))) ;
                    SPEX_CHECK(SPEX_mpz_mul(SPEX_1D(A, nnz, mpz),
                        SPEX_1D(A, nnz, mpz), SPEX_MPQ_NUM(B->v[j]->scale))) ;
                }
                nnz++;
            }
            A->p[j+1] = nnz;
        }
    }
    else
    {
        // if B is stored in row-wise format, find the column-wise nnz pattern
        // and # of nnz in B first
        rowcount  = (int64_t*) SPEX_calloc(B->m, sizeof(int64_t));
        if (!rowcount)
        {
            SPEX_FREE_ALL;
            return SPEX_OUT_OF_MEMORY;
        }

        // find the # of nz in each col of B
        for (i = 0; i < B->n; i++)
        {
            nnz += B->v[i]->nz;
            for (p = 0 ; p < B->v[i]->nz; p++)
            {
                j = B->v[i]->i[p];
                rowcount [j]++ ;
            }
        }
        // allocate space for A
        SPEX_CHECK(SPEX_matrix_allocate(&A, SPEX_CSC, SPEX_MPZ, B->m, B->n, nnz,
            false, true, option));

        // compute cumulative sum of rowcount to get the col pointer
        SPEX_CHECK(SPEX_cumsum(A->p, rowcount, B->m));

        // construct A from B
        for (i = 0 ; i < B->n ; i++)
        {
            SPEX_CHECK(SPEX_mpq_cmp_ui(&sgn, B->v[i]->scale, 1, 1));
            for (p = 0 ; p < B->v[i]->nz ; p++)
            {
                j = B->v[i]->i[p];
                int64_t rj = rowcount[j];
                A->i[ rj ] = (perm_inv == NULL) ? i : perm_inv[i] ;
                // A->x[ rowcount[j] ] = B->v[i]->x[p]*scale
                if (sgn == 0) // scale == 1
                {
                    SPEX_CHECK(SPEX_mpz_set(SPEX_1D(A, rj, mpz),
                        B->v[i]->x[p])) ;
                }
                else
                {
                    SPEX_CHECK(SPEX_mpz_divexact(SPEX_1D(A, rj, mpz),
                        B->v[i]->x[p], SPEX_MPQ_DEN(B->v[i]->scale))) ;
                    SPEX_CHECK(SPEX_mpz_mul(SPEX_1D(A, rj, mpz),
                        SPEX_1D(A, rj, mpz), SPEX_MPQ_NUM(B->v[i]->scale))) ;
                }
                rowcount[j] ++;
            }
        }
    }
    SPEX_CHECK(SPEX_mpq_set(A->scale, B->scale));

    (*A_handle) = A;
    SPEX_FREE_WORK;
    return SPEX_OK;
}
