//------------------------------------------------------------------------------
// SPEX_Update/SPEX_CSC_to_mat.c: convert a SPEX_matrix in CSC format to
// a SPEX_mat matrix that can be used in the SPEX LU update process
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// SPEX_CSC_to_mat create a SPEX_matr A from the given SPEX_matrix B which is
// in Compressed Sparse Column (CSC) format with entries in mpz_t type.  The
// entries of A will be in mpz_t type but the entries would be stored in either
// column form or row form, which is specified by A_Is_ColWise. Since the SPEX
// LU update process requires A = LDU, while the SPEX LU factorization gives
// PAQ = LDU, additional permutation will be applied.
//
// When converting A, we should use SPEX_CSC_to_mat(A, NULL, true, A1, option)
// and the resulting A will be in the same column and row order as A1.
//
// When converting L, we should use SPEX_CSC_to_mat(L, P, true, L1, option)
// and the resulting L will have same column permutation as L1 but same row
// order as A.
//
// When converting U, we should use SPEX_CSC_to_mat(U, Q, false, U1, option)
// and the resulting U will have same row permutation as U1, but same column
// order as A.

#define SPEX_FREE_WORK               \
    SPEX_FREE(rowcount);

#define SPEX_FREE_ALL                \
    SPEX_FREE_WORK;                  \
    SPEX_mat_free(&A);

#include "spex_update_internal.h"

SPEX_info SPEX_CSC_to_mat
(
    SPEX_mat **A_handle,          // converted SPEX_mat matrix
    const int64_t *perm,          // vector for permutation matrix, consider as
                                  // identity matrix if input as NULL
    const bool A_Is_ColWise,      // true if A->v[i] is the i-th col of A.
                                  // Otherwise, A->v[i] is the i-th row of A
    const SPEX_matrix *B,         // original matrix
    const SPEX_options *option
)
{
    SPEX_REQUIRE (B, SPEX_CSC, SPEX_MPZ) ;
    if (A_handle == NULL)   {return SPEX_INCORRECT_INPUT;}

    SPEX_info info;
    int64_t *rowcount = NULL;
    (*A_handle) = NULL;
    SPEX_mat *A = NULL;
    int64_t *Bp = B->p;

    int64_t i, j, p, Ap = 0;

    if (A_Is_ColWise) // if A will be stored in a column-wise format
    {
        // allocate space for A
        SPEX_CHECK(SPEX_mat_alloc(&A, B->n, B->m, true));

        for (j = 0; j < B->n; j++)
        {
            // reallocate space for each column of A
            SPEX_CHECK(SPEX_vector_realloc(A->v[j], Bp[j+1]-Bp[j]));
            Ap = 0;
            for (p = Bp[j]; p < Bp[j+1]; p++)
            {
                i = B->i[p];
                A->v[j]->i[Ap] = (perm == NULL) ? i : perm[i];
                // A->v[j]->x[Ap] = B->x[p]
                SPEX_CHECK(SPEX_mpz_set(A->v[j]->x[Ap], SPEX_1D(B, p, mpz)));
                Ap++;
            }
            A->v[j]->nz = Ap;
        }
    }
    else           // if A will be stored in row-wise format
    {
        // allocate space for A
        SPEX_CHECK(SPEX_mat_alloc(&A, B->m, B->n, true));

        // find the # of nz in each row of B
        rowcount  = (int64_t*) SPEX_calloc(B->m, sizeof(int64_t));
        if (!rowcount)
        {
            SPEX_FREE_ALL;
            return SPEX_OUT_OF_MEMORY;
        }
        for (p = 0; p < Bp[B->n]; p++)
        {
            rowcount[ B->i[p] ]++ ;
        }

        // reallocate space for each column of A
        for (i = 0; i < B->m; i++)
        {
            SPEX_CHECK(SPEX_vector_realloc(A->v[i], rowcount[i]));
        }

        // construct A from B
        for (j = 0; j < B->n; j++)
        {
            for (p = Bp[j]; p < Bp[j+1]; p++)
            {
                i = B->i[p];
                Ap = A->v[i]->nz++;
                A->v[i]->i[Ap] = (perm == NULL) ? j : perm[j] ;
                // A->v[i]->x[Ap] = B->x[p]
                SPEX_CHECK(SPEX_mpz_set(A->v[i]->x[Ap], SPEX_1D(B, p, mpz)));
            }
        }
    }
    SPEX_CHECK(SPEX_mpq_set(A->scale, B->scale));

    (*A_handle) = A;
    SPEX_FREE_WORK;
    return SPEX_OK;
}
