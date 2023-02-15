//------------------------------------------------------------------------------
// SPEX_Utilities/SPEX_factorization_convert: convert between updatable
// and non-updatable factorization.
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2020-2023, Christopher Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// The L and U factorization from SPEX_LU or the L from SPEX_Cholesky are
// all in SPEX_CSC format, and their columns and rows are permuted to be the
// same as the permuted matrix A(P,Q), and thus A(P,Q)=LD^(-1)U. However, all
// Update functions requires A=LD^(-1)U with L and/or U in SPEX_DYNAMIC_CSC
// format. This is the function to perform the in-place coversion for L and U
// so that it meet the requirement for Update function and vice versa.  If
// updatable == false, the returned factorization will be non-updatable with
// MPZ entries and CSC kind. Otherwise (updatable == true), the returned
// factorization will be updatable with MPZ entries and DYNAMIC_CSC kind.
//
// To help understand the process of the conversion, the steps to obtain a deep
// copy of converted factorization is provided below, which is slightly
// different from what is done in this function to obtain an in-place
// conversion:
//
// To get updatable (dynamic_CSC MPZ) from non-updatable (CSC MPZ) matrix:
//     1. performing SPEX_matrix_copy to get the L and/or UT (needs additional
//     transpose before copy) SPEX_DYNAMIC_CSC form.
//     2. permute row indices of L or UT such that A=LD^(-1)U, or equivalently
//     A_out->v[j]->i[p] = perm[A_in->v[j]->i[p]], where A is either L or UT.
//     3. canonicalize a SPEX_DYNAMIC_CSC matrix such that each column of
//     the input matrix have corresponding pivot as the first entry.
//
// To get non-updatable (CSC MPZ) from updatable (dynamic_CSC MPZ) matrix:
//     1. un-permute L and/or U by performing the inverse of the permutation.
//     2. call SPEX_matrix_copy to obtain L and/or U in the SPEX_CSC format, and
//     call SPEX_transpose for U if exists
//
// NOTE: if F->updatable == false upon input, F->L (and F->U if exists) must be
// CSC MPZ matrix, otherwise, SPEX_INCORRECT_INPUT will be returned. Likewise,
// if F->updatable == true upon input, F->L (and F->U if exists) must be
// dynamic_CSC MPZ matrix. In addition, both F->L and F->U (if exists) must not
// be shallow matrices. All SPEX functions output F with matrices in either of
// these two formats and non-shallow. Therefore, these input requirements can
// be met easily if users do not try to modify any individual component of F.
// The conversion is done in place.  In case of any error, the returned
// factorization should be considered as undefined.

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
            SPEX_vector_free(&(Mv[i]), option); \
        }                            \
        SPEX_FREE(Mv);               \
    }                                \
}

#include "spex_util_internal.h"

SPEX_info SPEX_factorization_convert
(
    SPEX_factorization F,       // The factorization to be converted
    bool updatable,             // if true: make F updatable
                                // if false: make non-updatable
    const SPEX_options option   // Command options
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------

    SPEX_info info = spex_factorization_basic_check (F);
    if (info != SPEX_OK) return (info);

    //--------------------------------------------------------------------------
    // quick return
    //--------------------------------------------------------------------------

    if (updatable == F->updatable)
    {
        // nothing to do
        return SPEX_OK ;
    }

    //--------------------------------------------------------------------------
    // initialize workspace
    //--------------------------------------------------------------------------

    SPEX_matrix L = F->L, U = F->U;
    int64_t i, j, p, m = L->m, n = L->n, mp, nnz = 0;
    SPEX_vector *Mv = NULL;
    mpz_t *Mx = NULL;
    int64_t *rowcount = NULL, *Mp = NULL, *Mi = NULL;

    // update the updatable flag
    F->updatable = updatable;

    //--------------------------------------------------------------------------
    // obtain/update Qinv_perm for updatable LU factorization
    //--------------------------------------------------------------------------
    if (F->kind == SPEX_LU_FACTORIZATION && updatable)
    {
        // Although Qinv_perm is NULL when F is created by factorizing matrix,
        // where F is initially not updatable, Qinv_perm is then created when
        // converted to updatable format and not deleted even when converted
        // back to non-updatable from updatable format.

        // Allocate space for Qinv_perm if it was NULL. Otherwise, Qinv_perm
        // is re-usable since the size is never changed.
        if (!(F->Qinv_perm))
        {
            F->Qinv_perm = (int64_t*) SPEX_malloc (n * sizeof(int64_t));
            if (!(F->Qinv_perm))
            {
                return SPEX_OUT_OF_MEMORY;
            }
        }
        for (i = 0; i < n; i++)
        {
            F->Qinv_perm[F->Q_perm[i]] = i;
        }
    }

    //--------------------------------------------------------------------------
    // obtain the desired format of L and/or U
    //--------------------------------------------------------------------------

    // The following converts either F->L or F->U between CSC MPZ
    // matrix and DYNAMIC_CSC MPZ matrix.

    // if updatable == true: This function converts F->L from CSC MPZ matrix to
    // Dynamic_CSC MPZ matrix in an updatable format, which means that rows of
    // L will be properly permuted with F->P_perm and the first entry of each
    // column of L will be the corresponding diagonal.  If F is an LU
    // factorization, then F->U will be also converted from CSC MPZ matrix to
    // Dynamic_CSC MPZ matrix in an updatable format, which means that U will
    // be firstly transposed to UT, and then the rows of UT will be properly
    // permuted  with F->Q_perm, and the first entry of each column of UT will
    // be the corresponding diagonal.

    // if updatable == false: This function converts F->L from Dynamic_CSC MPZ
    // matrix to CSC MPZ matrix in a non-updatable format, which means that the
    // permutation of the rows of L will be reset with F->Pinv_perm.  If F is
    // an LU factorization, F->U is also converted from Dynamic_CSC MPZ matrix
    // to CSC MPZ matrix in a non-updatable format, which means that U will be
    // firstly transposed to UT, and then the permutation of the rows of UT
    // will be reset with F->Qinv_perm.

    int64_t *perm;

    if (updatable)
    {

        //----------------------------------------------------------------------
        // convert a non-updatable LU or Cholesky factorization to updatable
        //----------------------------------------------------------------------

        //----------------------------------------------------------------------
        // allocate an array of n vetors to construct the DYNAMIC_CSC MPZ matrix
        //----------------------------------------------------------------------
        Mv = (SPEX_vector*) SPEX_calloc(n, sizeof(SPEX_vector));
        if (!Mv)
        {
            SPEX_FREE_ALL;
            return SPEX_OUT_OF_MEMORY;
        }

        for (i = 0; i < n; i++)
        {
            SPEX_CHECK(SPEX_vector_allocate(&(Mv[i]), 0, option));
        }

        //------------------------------------------------------------------
        // convert L from CSC MPZ to DYNAMIC_CSC MPZ
        //------------------------------------------------------------------

        int64_t *Lp = L->p;
        perm = F->P_perm;
        for (j = 0; j < n; j++)
        {
            // reallocate space for each column vector Mv[j]
            SPEX_CHECK(SPEX_vector_realloc(Mv[j], Lp[j+1]-Lp[j], option));
            mp = 0;
            bool found_diag = false;
            for (p = Lp[j]; p < Lp[j+1]; p++)
            {
                i = L->i[p];
                // permute entry indices
                Mv[j]->i[mp] = perm[i];
                // Mv[j]->x[mp] = L->x[p]
                SPEX_MPZ_SWAP(Mv[j]->x[mp], L->x.mpz[p]);

                // make sure the first entry is the diagonal
                if (!found_diag && i == j)
                {
                    found_diag = true;
                    if (mp != 0)
                    {
                        SPEX_MPZ_SWAP(Mv[j]->x[0],
                                                 Mv[j]->x[mp]);
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

        //----------------------------------------------------------------------
        // convert U if F is an LU factorization
        //----------------------------------------------------------------------

        if (F->kind == SPEX_LU_FACTORIZATION)
        {
            int64_t *Up = U->p;
            perm = F->Q_perm;

            //------------------------------------------------------------------
            // allocate an array of n vetors to construct the DYNAMIC_CSC MPZ
            // matrix
            //------------------------------------------------------------------
            Mv = (SPEX_vector*) SPEX_calloc(n, sizeof(SPEX_vector));
            if (!Mv)
            {
                SPEX_FREE_ALL;
                return SPEX_OUT_OF_MEMORY;
            }

            for (i = 0; i < n; i++)
            {
                SPEX_CHECK(SPEX_vector_allocate(&(Mv[i]), 0, option));
            }

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
                for (p = Up[j]; p < Up[j+1]; p++)
                {
                    i = U->i[p];
                    mp = Mv[i]->nz++;
                    Mv[i]->i[mp] = perm[j] ;
                    // Mv[i]->x[mp] = U->x[p]
                    SPEX_MPZ_SWAP(Mv[i]->x[mp], U->x.mpz[p]);

                    // Given U in CSC, getting UT in CSC is equivalent to
                    // obtaining U in CSR. Therefore, the pivot entry of j-th
                    // column of U (in CSC) is always the first entry
                    // being inserted to the j-th row of U (in CSR).
                    // In fact, all entries in each vector will be in order.
                    // Thus, i == j is true if and only if mp == 0 is true.
                    ASSERT ((i == j) == (mp == 0)) ;
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
    else
    {

        //----------------------------------------------------------------------
        // convert a updatable LU or Cholesky factorization to non-updatable
        //----------------------------------------------------------------------

        int sgn;

        //------------------------------------------------------------------
        // convert L from DYNAMIC_CSC MPZ to CSC MPZ
        //------------------------------------------------------------------

        perm = F->Pinv_perm;
        // compute number of nonzero in L
        for (j = 0; j < n; j++)
        {
            nnz += L->v[j]->nz;
        }

        // allocate space for Mi, Mp and Mx
        Mp = (int64_t *) SPEX_calloc (n+1, sizeof (int64_t));
        Mi = (int64_t *) SPEX_calloc (nnz, sizeof (int64_t));
        Mx = spex_create_mpz_array (nnz);
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
            SPEX_MPQ_CMP_UI(&sgn, L->v[j]->scale, 1, 1);
            for (p = 0 ; p < L->v[j]->nz ; p++)
            {
                i = L->v[j]->i[p];
                Mi[mp] = perm[i] ;
                if (sgn != 0) // scale != 1
                {
                    // apply scale to L->v[j]->x[p]
                    SPEX_MPZ_DIVEXACT(L->v[j]->x[p],
                        L->v[j]->x[p], SPEX_MPQ_DEN(L->v[j]->scale));
                    SPEX_MPZ_MUL(L->v[j]->x[p],
                        L->v[j]->x[p], SPEX_MPQ_NUM(L->v[j]->scale));
                }
                SPEX_MPZ_SWAP(Mx[mp], L->v[j]->x[p]);
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

        //----------------------------------------------------------------------
        // convert U if F is an LU factorization
        //----------------------------------------------------------------------

        if (F->kind == SPEX_LU_FACTORIZATION)
        {
            perm = F->Qinv_perm;
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
                    j = perm[U->v[i]->i[p]];
                    rowcount [j]++ ;
                }
            }

            // allocate space for Mi, Mp and Mx
            Mp = (int64_t *) SPEX_calloc (n+1, sizeof (int64_t));
            Mi = (int64_t *) SPEX_calloc (nnz, sizeof (int64_t));
            Mx = spex_create_mpz_array (nnz);
            if (!Mp || !Mi || !Mx)
            {
                SPEX_FREE_ALL;
                return SPEX_OUT_OF_MEMORY;
            }

            // compute cumulative sum of rowcount to get the col pointer
            SPEX_CHECK(spex_cumsum(Mp, rowcount, m));

            // construct Mi, Mp and Mx from U
            for (i = 0 ; i < n ; i++)
            {
                SPEX_MPQ_CMP_UI(&sgn, U->v[i]->scale, 1, 1);
                for (p = 0 ; p < U->v[i]->nz ; p++)
                {
                    j = perm[U->v[i]->i[p]];
                    mp = rowcount[j] ++;
                    Mi[ mp ] = i ;
                    if (sgn != 0) // scale != 1
                    {
                        // apply scale to U->v[i]->x[p]
                        SPEX_MPZ_DIVEXACT(U->v[i]->x[p],
                            U->v[i]->x[p], SPEX_MPQ_DEN(U->v[i]->scale));
                        SPEX_MPZ_MUL(U->v[i]->x[p],
                            U->v[i]->x[p], SPEX_MPQ_NUM(U->v[i]->scale));
                    }
                    SPEX_MPZ_SWAP(Mx[mp], U->v[i]->x[p]);
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

