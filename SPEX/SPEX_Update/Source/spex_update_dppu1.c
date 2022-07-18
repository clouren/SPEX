//------------------------------------------------------------------------------
// SPEX_Update/spex_update_dppu1.c: perform diagonal permutation pivot update
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Purpose: This function is called to perform diagonal permutation pivot update
// when the submatrix (formed by rows and columns k to ks) has the following
// pattern
//       x 0 0 0 0       <- row k
//       0 x . . 0
//       0 . x . 0
//       0 . . x 0
//       . 0 0 0 x       <- row ks
// The caller of this function always finds the biggest index from k+1:n-1 as ks
// such that the submatrix has the above pattern. The patterns of ks+1 row and
// column are unkown but must have at least one nonzero in U(k:ks,ks+1) or
// L(ks+1, k+1:ks).
//
// This function will swap rows and columns k and ks in L and U. Noted that the
// rows of L and columns of U are permuted implicitly via the permutation
// matrices based on P and Q.

#define SPEX_FREE_ALL                \
{                                    \
    SPEX_MPZ_CLEAR(Lksk);            \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPZ_CLEAR(tmpz);            \
}

#include "spex_update_internal.h"

#define SL(k) (L->v[(k)]->scale)
#define SU(k) (U->v[(k)]->scale)

SPEX_info spex_update_dppu1
(
    SPEX_matrix *L,     // matrix L
    SPEX_matrix *U,     // matrix U
    SPEX_matrix *rhos,// array of scaled pivots
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    int64_t *inext,  // the index of first off-diag entry in col k of L
    int64_t *h,      // allocated vector that can be used for history vector.
                     // All entries are maintained to be >= -1
    int64_t *Q,      // column permutation
    int64_t *Q_inv,  // inverse of column permutation
    int64_t *P,      // row permutation
    int64_t *P_inv,  // inverse of row permutation
    const int64_t k,   // current column index 0 <= k < n
    const int64_t ks,  // index of the diagonal to be swapped with, [0,n)
    const SPEX_options *option
)
{
    // initialize workspace
    SPEX_info info;
    int sgn, Lksk_sgn;
    int64_t pk, ck, pks, cks, j, n = U->n;
    int64_t Qk = Q[k], Pk = P[k], Qks, Pks;
    mpz_t *sd = rhos->x.mpz;
    SPEX_vector *v;

    mpq_t pending_scale;
    SPEX_MPQ_SET_NULL(pending_scale);
    mpz_t Lksk, tmpz; SPEX_MPZ_SET_NULL(Lksk); SPEX_MPZ_SET_NULL(tmpz);
    SPEX_CHECK(SPEX_mpq_init(pending_scale));
    SPEX_CHECK(SPEX_mpz_init(Lksk));
    SPEX_CHECK(SPEX_mpz_init(tmpz));

    // -------------------------------------------------------------------------
    // handle the special case when swapping with the inserted column. Since it
    // is only in the k-th IPGE iteration, there is no need to perform
    // backtracking for the inserted column. Therefore, only need to perform
    // backtracking for U(n-1, Q(n-1)).
    // -------------------------------------------------------------------------
    if (ks == n)
    {
        Qks = Q[n-1];
        Pks = P[n-1];

#ifdef SPEX_DEBUG
        // U(k,Q(ks)) must be nnz, otherwise it would cause singularity, since
        // U(k,Q(k+1:n)) will be all zeros
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row->x[Qks]));
        ASSERT(sgn != 0);
#endif

        // ---------------------------------------------------------------------
        // backtrack U(n-1,Q(n-1)) which is initially sd[n-1] with pending
        // scalors applied
        // ---------------------------------------------------------------------
        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpz_mul(U->v[n-1]->x[0], sd[n-1], sd[k-1]));
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_set(U->v[n-1]->x[0], sd[n-1]));
        }

        // L(P[ks],k) should be found at L(P[n-1],k) instead
        SPEX_CHECK(SPEX_mpz_sgn(&(Lksk_sgn), Lk_dense_col->x[Pks]));
#ifdef SPEX_DEBUG
        printf("using dppu1 swapping k(%ld) and n(%ld), L(ks,k)%s=0\n",k,ks,Lksk_sgn==0?"=":"!");
#endif
        if (Lksk_sgn == 0)
        {
            SPEX_CHECK(SPEX_mpz_divexact(U->v[n-1]->x[0],
                                         U->v[n-1]->x[0], sd[n-2]));
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_fdiv_q(U->v[n-1]->x[0],
                                       U->v[n-1]->x[0], sd[n-2]));
            // L(P(n-1),k) = L(P(n-1),k)*S(1,k), which should be integer
            SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[Pks],
                       Lk_dense_col->x[Pks], SPEX_MPQ_DEN(SL(k))));
            SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[Pks],
                       Lk_dense_col->x[Pks], SPEX_MPQ_NUM(SL(k))));

            // tmpz = ceil(U(k,Q(n-1))*L(P(n-1),k)/U(k,Q(k))
            SPEX_CHECK(SPEX_mpz_mul(tmpz, Uk_dense_row->x[Qks],
                                    Lk_dense_col->x[Pks]));
            SPEX_CHECK(SPEX_mpz_cdiv_q(tmpz, tmpz, Uk_dense_row->x[Qk]));

            // U(n-1,Q(n-1)) = U(n-1,Q(n-1))+tmpz
            SPEX_CHECK(SPEX_mpz_add(U->v[n-1]->x[0], U->v[n-1]->x[0], tmpz));

            // reset Lk_dense_col->x[P[n-1]]=0
            SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col->x[Pks], 0));
        }
        // S(2,n-1) is considered as 1 beyond this point, even though it is
        // not explicit set to 1 here. It will be updated after calling this
        // function.

        // ---------------------------------------------------------------------
        // scale entries in frames k+1:n-2
        // ---------------------------------------------------------------------
        // since the value in Uk_dense_row[Q[k]] will not be used, we use it to
        // hold the original value of sd[k] before swapping columns and rows of
        // k and n-1. Then we set sd[k] to the new pivot of column k of L
        // vk[P[k]], which is keep as first entry in the nnz list.
        SPEX_CHECK(SPEX_mpz_swap(Uk_dense_row->x[Qk], sd[k]));
        SPEX_CHECK(SPEX_mpz_set(sd[k], L->v[k]->x[0]));

        if (n > k+2) // n-1 > k+1
        {
            // pending_scale = sd(k)/Uk_dense_row[Q[k]]
            SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k]));
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, Uk_dense_row->x[Qk]));
            SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

            for (j = k+1; j < n-1; j++)
            {
                // S(:,k+1:n-2) = S(:,k+1:n-2)*pending_scale;
                SPEX_CHECK(SPEX_mpq_mul(SL(j), SL(j), pending_scale));
                SPEX_CHECK(SPEX_mpq_mul(SU(j), SU(j), pending_scale));
                // sd(k+1:n-2) = sd(k+1:n-2)*pending_scale;
                SPEX_CHECK(SPEX_mpz_divexact(sd[j],
                                        sd[j], SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(sd[j],
                                        sd[j], SPEX_MPQ_NUM(pending_scale)));
            }
        }

        // ---------------------------------------------------------------------
        // Columns of k and n will be swapped after calling this function, we
        // only need to swap rows of k and n-1
        // ---------------------------------------------------------------------
        // swap rows k and n-1 of U           % O(1) time
        v = U->v[k];       U->v[k] = U->v[n-1];    U->v[n-1] = v;

        // update row permutation to swap rows of L implicitly
        P[k] = Pks;          P[n-1] = Pk;
        P_inv[Pks] = k;   P_inv[Pk] = n-1;

        // ---------------------------------------------------------------------
        // In order to get sd(n-1), we just need to apply skip scale for
        // U(k, Q[n-1]) and perform IPGE update. Then use it to update
        // L(P[n-1],n-1), U(n-1,Q[n-1]) and S for frame n-1
        // ---------------------------------------------------------------------
        // get the scale for IPGE update: pending_scale = sd(n-2)/sd(k-1);
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[n-2]));
        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k-1]));
            SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        }
        // pending_scale *= S(2,n-1), which was S(2,k) before swapping
        SPEX_CHECK(SPEX_mpq_mul(pending_scale, pending_scale, SU(n-1)));

        // sd[n-1] = U(k, Q[n-1]) * pending_scale
        SPEX_CHECK(SPEX_mpz_swap(sd[n-1], Uk_dense_row->x[Qks]));
        SPEX_CHECK(SPEX_mpz_divexact(sd[n-1], sd[n-1],
                                     SPEX_MPQ_DEN(pending_scale)));
        SPEX_CHECK(SPEX_mpz_mul(sd[n-1], sd[n-1],
                                     SPEX_MPQ_NUM(pending_scale)));

        // set U(n-1,Q(n-1))=L(P(n-1),n-1)=sd[n-1]
        SPEX_CHECK(SPEX_mpz_set (L->v[n-1]->x[0], sd[n-1]));
        SPEX_CHECK(SPEX_mpz_set (U->v[n-1]->x[0], sd[n-1]));
        U->v[n-1]->i[0] = Qks;
        U->v[n-1]->nz = 1;
        L->v[n-1]->i[0] = Pk;
        L->v[n-1]->nz = 1;

        // S(:,n-1) = [1;1]
        SPEX_CHECK(SPEX_mpq_set_ui(SL(n-1), 1, 1));
        SPEX_CHECK(SPEX_mpq_set_ui(SU(n-1), 1, 1));

        // reset nnz entries in Uk_dense_row to 0
        SPEX_CHECK(SPEX_mpz_set_ui(Uk_dense_row->x[Qk], 0));
        SPEX_CHECK(SPEX_mpz_set_ui(Uk_dense_row->x[Qks], 0));

        SPEX_FREE_ALL;
        return SPEX_OK;
    }

    Qks = Q[ks];
    Pks = P[ks];
    // -------------------------------------------------------------------------
    // perform backtracking for row ks of U
    // -------------------------------------------------------------------------
    // find the scale for backtracking: pending_scale = sd(k-1)/sd(ks-1)
    if (k == 0)
    {
        SPEX_CHECK(SPEX_mpq_set_ui(pending_scale, 1, 1));
    }
    else
    {
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k-1]));
    }
    SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[ks-1]));
    // remove common factor in mpq_den and mpq_num
    SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

    // S(:, ks) = pending_scale*S(:,ks)
    SPEX_CHECK(SPEX_mpq_mul(SL(ks), SL(ks), pending_scale));
    SPEX_CHECK(SPEX_mpq_mul(SU(ks), SU(ks), pending_scale));

    SPEX_CHECK(SPEX_mpz_sgn(&(Lksk_sgn), Lk_dense_col->x[Pks]));
#ifdef SPEX_DEBUG
    printf("using dppu1 swapping k(%ld) and ks(%ld), L(ks,k)%s=0\n",
        k, ks, Lksk_sgn==0?"=":"!");
#endif
    if (Lksk_sgn == 0) // L(P(ks),k) == 0
    {
        ASSERT (*inext != ks);
        // sd(ks) = sd(ks)*pending_scale
        SPEX_CHECK(SPEX_mpz_divexact(sd[ks],
                                     sd[ks], SPEX_MPQ_DEN(pending_scale)));
        SPEX_CHECK(SPEX_mpz_mul(sd[ks], sd[ks], SPEX_MPQ_NUM(pending_scale)));
    }
    else
    {
        ASSERT (*inext == ks);
        // assign Lksk = L(P(ks),k)*S(1,k), which should be integer
        SPEX_CHECK(SPEX_mpz_divexact(Lksk, Lk_dense_col->x[Pks],
                                     SPEX_MPQ_DEN(SL(k))));
        SPEX_CHECK(SPEX_mpz_mul(Lksk, Lksk, SPEX_MPQ_NUM(SL(k))));

        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        // backtracking jumbled sparse row ks of U using scattered row k of U.
        // Assuming explicit zeros in U(ks, :) resulted from exact cancellation
        // in IPGE update were not removed (SLIP LU keeps those zeros in output
        // L and U), nonzero pattern of U(k,Q(k+1:n+1)) should be a subset of
        // U(ks,:).
        // However, we remove every explicit zero we found for both IPGE and
        // backtracking. We will need to insert new entry to row ks after
        // backtracking. To this end, we need to iterate across all nonzeros
        // in Uk_dense_row to find if any column index of nonzero is untouched,
        // then a new nonzero should be added.
        // U(ks,cks) = U(k,cks)*Lksk/U(k,Q(k))
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        int64_t count = 0; // # of nnz in row k of U, explicit zeros excluded
        for (pk = 0; pk < Uk_dense_row->nz; pk++)
        {
            ck = Uk_dense_row->i[pk];
            SPEX_CHECK(SPEX_mpz_sgn (&sgn, Uk_dense_row->x[ck]));
            if (sgn != 0)
            {
                h[ck] = -2;
                count ++;
            }
        }
        // exclude diagonal
        h[Qk] = -1;
        count --;

        // start backtracking
        pks = 0;
        while (pks < U->v[ks]->nz)
        {
            // column index in row ks of U
            cks = U->v[ks]->i[pks];

            // U(ks,cks) = U(ks,cks)*S(2,ks)+U(k,cks)*Lksk/U(k,Q(k))
            if (h[cks] == -2)  // U(k,cks) != 0
            {
                // tmpz = ceil(U(k,cks)*Lksk/U(k,Q(k))
                SPEX_CHECK(SPEX_mpz_mul(tmpz, Uk_dense_row->x[cks], Lksk));
                SPEX_CHECK(SPEX_mpz_cdiv_q(tmpz, tmpz, Uk_dense_row->x[Qk]));

                // U(ks,cks) = floor(U(ks,cks)*S(2,ks))
                SPEX_CHECK(SPEX_mpz_mul(U->v[ks]->x[pks], U->v[ks]->x[pks],
                                        SPEX_MPQ_NUM(SU(ks))));
                SPEX_CHECK(SPEX_mpz_fdiv_q(U->v[ks]->x[pks], U->v[ks]->x[pks],
                                        SPEX_MPQ_DEN(SU(ks))));

                // U(ks,cks) = U(ks,cks)+tmpz
                SPEX_CHECK(SPEX_mpz_add(U->v[ks]->x[pks],
                                        U->v[ks]->x[pks], tmpz));

                h[cks] = -1;
                count--;

                // remove this entry if it becomes zero
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[ks]->x[pks]));
                if (sgn == 0)
                {
                    U->v[ks]->nz --;
                    U->v[ks]->i[pks] = U->v[ks]->i[U->v[ks]->nz];
                    SPEX_CHECK(SPEX_mpz_swap(U->v[ks]->x[pks],
                                             U->v[ks]->x[U->v[ks]->nz]));
                    continue;
                }
            }
            else           // U(k,cks) == 0
            {
                // U(ks,cks) = U(ks,cks)*S(2,ks)
                SPEX_CHECK(SPEX_mpz_divexact(U->v[ks]->x[pks], U->v[ks]->x[pks],
                                        SPEX_MPQ_DEN(SU(ks))));
                SPEX_CHECK(SPEX_mpz_mul(U->v[ks]->x[pks], U->v[ks]->x[pks],
                                        SPEX_MPQ_NUM(SU(ks))));

                // update sd[ks] = U(ks, Q(ks))
                if (cks == Qks)
                {
                    SPEX_CHECK(SPEX_mpz_set(sd[ks], U->v[ks]->x[pks]));
                }
            }
            pks ++;
        }
        // continue backtracking for potential fillin
        if (count > 0)
        {
            // allocate additional space if needed
            if (U->v[ks]->nz+count > U->v[ks]->nzmax)
            {
                SPEX_CHECK(SPEX_vector_realloc(U->v[ks], U->v[ks]->nz+count,
                    option));
            }
            pks = U->v[ks]->nz;
            for (pk = 0; count > 0 && pk < Uk_dense_row->nz; pk++)
            {
                ck = Uk_dense_row->i[pk];
                if (h[ck] == -2)
                {
                    U->v[ks]->i[pks] = ck;
                    // U(ks,ck) = U(k,ck)*Lksk/U(k,Q(k))
                    SPEX_CHECK(SPEX_mpz_mul(U->v[ks]->x[pks],
                        Uk_dense_row->x[ck], Lksk));
                    SPEX_CHECK(SPEX_mpz_divexact(U->v[ks]->x[pks],
                        U->v[ks]->x[pks], Uk_dense_row->x[Qk]));
                    pks++;

                    h[ck] = -1;
                    count--;
                }
            }
            U->v[ks]->nz = pks;
        }

        // reset S(2,ks)
        SPEX_CHECK(SPEX_mpq_set_ui(SU(ks), 1, 1));

        // Mathematically, we should insert new entry at U(ks, Q[k]) and swap
        // its value with Lksk. However, since the value of this entry will not
        // be used beyond this point, and column Q[k] of U will be deleted when
        // finished, we will skipped adding U(ks, Q[k]) here.
    }
    
    // ------------------------------------------------------------------------
    // swap rows and columns of k and ks
    // ------------------------------------------------------------------------
    // swap columns k and ks of L        % O(1) time
    v = L->v[k];       L->v[k] = L->v[ks];    L->v[ks] = v;
    // swap rows k and ks of U           % O(1) time
    v = U->v[k];       U->v[k] = U->v[ks];    U->v[ks] = v;

    // update row/column permutation to swap rows of L and cols of U implicitly
    Q[k] = Qks;           Q[ks] = Qk;
    Q_inv[Qks] = k;   Q_inv[Qk] = ks;
    P[k] = Pks;           P[ks] = Pk;
    P_inv[Pks] = k;   P_inv[Pk] = ks;

    // swap entries in sd
    SPEX_CHECK(SPEX_mpz_swap(sd[k], sd[ks]));

    // ------------------------------------------------------------------------
    // scale entries in frames k+1:ks-1
    // ------------------------------------------------------------------------
    if (ks > k+1)
    {
        // get the scale for entries between frames k and ks % O(1) time
        // pending_scale = sd(k)/sd (ks);
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k]));
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[ks]));
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

        for (j = k+1; j < ks; j++)
        {
            // S(:,k+1:ks-1) = S(:,k+1:ks-1)*pending_scale;
            SPEX_CHECK(SPEX_mpq_mul(SL(j), SL(j), pending_scale));
            SPEX_CHECK(SPEX_mpq_mul(SU(j), SU(j), pending_scale));
            // sd(k+1:ks-1) = sd(k+1:ks-1)*pending_scale;
            SPEX_CHECK(SPEX_mpz_divexact(sd[j],
                                    sd[j], SPEX_MPQ_DEN(pending_scale)));
            SPEX_CHECK(SPEX_mpz_mul(sd[j],
                                    sd[j], SPEX_MPQ_NUM(pending_scale)));
        }
    }

    // ------------------------------------------------------------------------
    // perform IPGE for frame ks, whose entries are in Lk_dense_col and
    // Uk_dense_row.
    // NOTE: in the last iteration, i.e., swapping column k with last column,
    // the IPGE update column k of L is useless, since the updated column will
    // be deleted and replaced. Therefore, its IPGE update in the last
    // iteration can be treated same as row k of U for better efficiency.
    // ------------------------------------------------------------------------
    // If this is the last iteration, then we don't need to perform the
    // remaining IPGE update, since L(:,k) will be deleted and updated with
    // inserted column.
    // However, the current heurestic used in SPEX_Update_LU_ColRep will always
    // handle such case with cppu instead. In case the heurestic is changed in
    // the furture, and thus this case become possible for dppu1 (the following
    // assert is triggered), enable the following if clause.
    ASSERT(ks != n-1); // This is used to detect if this case becomes possible.
    /*if (ks == n-1)
    {
        SPEX_FREE_ALL;
        return SPEX_OK;
    }*/

    // Since L(P[ks],k) will be 0 after swapping, the IPGE update for row k of
    // U can be done by multiplying with sd(ks-1)/sd(k-1).
    // get the scale for IPGE update: pending_scale = sd(ks-1)/sd (k-1);
    SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[ks-1]));
    if (k > 0)
    {
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k-1]));
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
    }
    // S(2, ks) = S(2, ks)*pending_scale
    SPEX_CHECK(SPEX_mpq_mul(SU(ks), SU(ks), pending_scale));

    // sd(ks) = sd(ks)*pending_scale;
    SPEX_CHECK(SPEX_mpz_divexact(sd[ks], sd[ks], SPEX_MPQ_DEN(pending_scale)));
    SPEX_CHECK(SPEX_mpz_mul(sd[ks], sd[ks], SPEX_MPQ_NUM(pending_scale)));

    if (Lksk_sgn == 0)
    {
#ifndef SPEX_DEBUG
        // disable this in the debug mode for simpler verification
        if (Lk_dense_col->nz == 1)
        {
            SPEX_FREE_ALL;
            return SPEX_OK;
        }
#endif
        // IPGE update for col k of L can be achieved by setting
        // S(1, ks) = S(1, ks)*pending_scale
        SPEX_CHECK(SPEX_mpq_mul(SL(ks), SL(ks), pending_scale));
    }
    else
    {
        Pk = Pks;  // Pk = P[k]
        // L(ck,ks) = (L(ck,ks)*L(P(k),k)-L(ck,k)*L(P(k),ks))/sd(k-1);
        // This IPGE euqation indicates that using L(P(k), ks) without applying
        // any pending scale factor to could keep it in the same pending scale
        // factor as the rest of entries in column k of L.  Therefore, the
        // skipped scaling for col k of L can still be skipped when performing
        // IPGE update. However, the result could be then non-integer. To avoid
        // that, we first set
        // S(1,ks) = S(1,ks)/sd[k-1]
        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpz_mul(SPEX_MPQ_DEN(SL(ks)),
                                    SPEX_MPQ_DEN(SL(ks)), sd[k-1]));
            SPEX_CHECK(SPEX_mpq_canonicalize(SL(ks)));
        }

        // In addition, we will skip applying pending scale to column k of L.
        // According to the IPGE update equation, we can see that the scale for
        // column ks of L S(1, ks) needs to multiply with S(1,k) after the IPGE
        // update.  We will do it in advance at this point before performing
        // the IPGE update, and divide all entries in column ks of L by the
        // denominator of S(1,ks), which will keep the result still in the
        // integer domain and also make it as small as possible.
        // S(1,ks) = S(1,ks)*S(1,k)
        SPEX_CHECK(SPEX_mpq_mul(SL(ks), SL(ks), SL(k)));

        // initialize history vector
        for (pks = 0; pks < Lk_dense_col->nz; pks++)
        {
            cks = Lk_dense_col->i[pks];
            // formally, we should set h[cks] = SPEX_FLIP(k-1), so we will know
            // the entries in L(:,ks) are in (k-1)-th IPGE iteration. However,
            // since we need to perform only one IPGE iteration, we just need
            // to know whether the corresponding entry is updated. Therefore,
            // the initialization for history vector is set as
            h[cks] = -2; // only entry in the nnz patter has h < -1
        }
        // # of nnz in column ks of L
        int64_t Lks_nz = Lk_dense_col->nz;

        // finish the IPGE update by performing
        // L(ck,ks) = L(ck,ks)*L(P(k),k)-L(ck,k)*L(P(k),ks))
        // NOTE: This will cause fillin in the ks(th) column of L,
        //       Since there is no subset relation between nnz pattern in
        //       L(:,ks) and L(:, k). Both could have explicit zero(s).
        //       L(:,k) can be jumbled.
        *inext = n;
        for (pk = 1 /*exclude L(P[k],k)*/; pk < L->v[k]->nz; pk++)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[k]->x[pk]));
            if (sgn == 0)  {continue;}// skip if L(ck, k) == 0

            // row index in column k of L
            ck = L->v[k]->i[pk];

            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[ck]));
            if (sgn != 0) // L(ck, ks) != 0
            {
                // L(ck, ks) = L(ck, ks) * L(P[k], k)
                SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ck],
                                        Lk_dense_col->x[ck],
                                        L->v[k]->x[0]));
            }
            else if (h[ck] >= -1) // this entry was not in nnz pattern
            {
                // insert new entry in the nonzero pattern
                ASSERT(Lks_nz <= Lk_dense_col->nzmax);
                Lk_dense_col->i[Lks_nz] = ck;
                Lks_nz++;
            }

            // L(ck, ks) -= L(P(k),ks)*L(ck,k) 
            SPEX_CHECK(SPEX_mpz_submul(Lk_dense_col->x[ck],
                                       Lk_dense_col->x[Pk], L->v[k]->x[pk]));

            // Check if resulting L(ck,ks) == 0
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[ck]));
            if (sgn != 0) // L(ck, ks) != 0
            {
                // only perform L(ck, ks) /= den(SL(ks)), and set
                // SL(ks) = num(SL(ks)) when finished. In this way, L(ck, ks)
                // will be smaller in size.
                SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ck],
                                             Lk_dense_col->x[ck],
                                             SPEX_MPQ_DEN(SL(ks))));

                // check if this will be the first off diagonal entry in L(P,ks)
                int64_t real_ck = P_inv[ck];
                if (real_ck < *inext)
                {
                    ASSERT(real_ck > ks);
                    // inext is the row index of the found first off-diagonal
                    // entry in L(P,ks)
                    *inext = real_ck;
                }
            }
            
            // update h[ck] to mark Lk_dense_col[ck] need no further update
            h[ck] = -1;
        }
        pks = 0;
        while (pks < Lks_nz)       // iterate across all nnz
        {
            // row index in column ks of L
            cks = Lk_dense_col->i[pks];

            if (h[cks] < -1) //only need to update entries that were not updated
            {
                // reset history vector to any value >= -1
                h[cks] = -1;
                // L(P(k), ks) should be removed from nnz pattern
                if (cks == Pk)
                {
                    // update the number of nnz
                    Lks_nz--;
                    // move the row index of last nonzero to current position
                    Lk_dense_col->i[pks] = Lk_dense_col->i[Lks_nz];
                    SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col->x[cks], 0));
                    continue;
                }

                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[cks]));
                if (sgn != 0) // L(cks, ks) != 0
                {
                    // L(cks,ks) = (L(cks,ks)*L(P(k),k))/den(S(1,ks));
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[cks],
                                            Lk_dense_col->x[cks],
                                            L->v[k]->x[0]));
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[cks],
                                            Lk_dense_col->x[cks],
                                            SPEX_MPQ_DEN(SL(ks))));
                    // check if this will be the 1st off-diag entry in L(P,ks)
                    int64_t real_cks = P_inv[cks];
                    if (real_cks < *inext && real_cks > ks)
                    {
                        // inext is the row index of the found first
                        // off-diagonal entry in L(P,ks)
                        *inext = real_cks;
                    }
                }
            }
            pks++;
        }
        // update the number of nnz
        Lk_dense_col->nz = Lks_nz;
        
        // skip the rest of IPGE iterations
        // pending_scale = sd(ks-1)/sd(k);
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[ks-1]));
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k]));
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        // set denominator of S(1,ks) = 1
        SPEX_CHECK(SPEX_mpz_set_ui(SPEX_MPQ_DEN(SL(ks)), 1));
        // S(2, ks) = S(2, ks)*pending_scale;
        SPEX_CHECK(SPEX_mpq_mul(SL(ks), SL(ks), pending_scale));
    }
    SPEX_FREE_ALL;
    return SPEX_OK;
}
