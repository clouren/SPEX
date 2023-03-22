//------------------------------------------------------------------------------
// SPEX_Update/spex_update_dppu2: perform diagonal permutation pivot update
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2023, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// Purpose: This function is called to perform diagonal permutation pivot update
// when the submatrix (formed by columns k to ks) has the following
// pattern
//       . . . . .  (row 1)
//       . . . . .  (rows 2 to k-1)
//       x . . . .  (row k)
//       0 x . . 0  (row k+1)
//       0 . x . 0  ( .... )
//       0 . . x 0  (row ks-1)
//       0 0 0 0 x  (row ks)
//       0 . . . .  (row ks+1 to n-1)
//       0 . . . .  (row n)
// This function will swap rows and columns k and ks in L and U. Noted that the
// rows of L and columns of U are permuted implicitly via the permutation
// matrices based on P and Q.

#define SPEX_FREE_ALL                \
{                                    \
    SPEX_MPQ_CLEAR(pending_scale);   \
}

#include "spex_update_internal.h"

#define SL(k) (L->v[(k)]->scale)
#define SU(k) (U->v[(k)]->scale)

SPEX_info spex_update_dppu2
(
    SPEX_matrix L,      // matrix L
    SPEX_matrix U,      // matrix U
    SPEX_matrix rhos,   // array of scaled pivots
    spex_scattered_vector Lk_dense_col, // scattered column k of L
    spex_scattered_vector Uk_dense_row, // scattered column k of U
    int64_t *jnext,     // the index of first off-diag entry in row k of U
    int64_t *h,         // allocated vector that can be used for history vector.
                        // All entries are maintained to be >= -1
    int64_t *Q,         // column permutation
    int64_t *Q_inv,     // inverse of column permutation
    int64_t *P,         // row permutation
    int64_t *P_inv,     // inverse of row permutation
    const int64_t k,    // current column index 0 <= k < n
    const int64_t ks,   // index of the diagonal to be swapped with, [0,n)
    const SPEX_options option
)
{

    // initialize workspace
    SPEX_info info;
    int sgn;
    int64_t j, n = U->n, tmp_ks = SPEX_MIN(ks, n-1);
    int64_t Qk = Q[k], Pk = P[k], Qks = Q[tmp_ks], Pks = P[tmp_ks], Qj;
    mpz_t *sd = rhos->x.mpz;
    SPEX_vector v; // used to switch vectors, no need to allocate nor free

    mpq_t pending_scale;
    SPEX_MPQ_SET_NULL(pending_scale);
    SPEX_MPQ_INIT(pending_scale);

#ifdef SPEX_DEBUG
    printf("using dppu2 with k(%ld) and ks(%ld)\n",k,ks);
#endif

    if (ks == n)
    {
        //----------------------------------------------------------------------
        // backtrack row n-1 of frame matrix, which only has 1 entry U(n-1,n-1)
        //----------------------------------------------------------------------
        if (k > 0)
        {
            SPEX_MPZ_MUL(U->v[n-1]->x[0], sd[n-1], sd[k-1]);
            SPEX_MPZ_DIVEXACT(U->v[n-1]->x[0],
                                         U->v[n-1]->x[0], sd[n-2]);
        }
        else
        {
            SPEX_MPZ_DIVEXACT(U->v[n-1]->x[0], sd[n-1], sd[n-2]);
        }
        // S(2,n-1) is considered as 1 beyond this point, even though it is
        // not explicit set to 1 here. It will be updated after calling this
        // function

        //----------------------------------------------------------------------
        // update entries in frames between k and n-1
        //----------------------------------------------------------------------
        // since the value in Uk_dense_row[Q[k]] will not be used, we use it to
        // hold the original value of sd[k] before swapping columns and rows of
        // k and n-1. Then we set sd[k] to the new pivot of column k of L
        // vk[P[k]], which is kept as first entry in the nnz list.
        SPEX_MPZ_SWAP(Uk_dense_row->x[Qk], sd[k]);
        SPEX_MPZ_SET(sd[k], L->v[k]->x[0]);

        if (n > k+2) // n-1 > k+1, i.e., row k is not the 2nd last row
        {
            // pending_scale = sd(k)/Uk_dense_row[Q[k]]
            SPEX_MPQ_SET_Z(pending_scale, sd[k]);
            SPEX_MPQ_SET_DEN(pending_scale, Uk_dense_row->x[Qk]);
            SPEX_MPQ_CANONICALIZE(pending_scale);

            // scale entries in frames k+1:n-2
            for (j = k+1; j < n-1; j++)
            {
                // S(:,j) = S(:,j)*pending_scale;
                SPEX_MPQ_MUL(SL(j), SL(j), pending_scale);
                SPEX_MPQ_MUL(SU(j), SU(j), pending_scale);
                // sd(j) = sd(j)*pending_scale;
                SPEX_MPZ_DIVEXACT(sd[j],
                                        sd[j], SPEX_MPQ_DEN(pending_scale));
                SPEX_MPZ_MUL(sd[j],
                                        sd[j], SPEX_MPQ_NUM(pending_scale));
            }
        }

        // reset Uk_dense_row->x[Q[k]]=0
        SPEX_MPZ_SET_UI(Uk_dense_row->x[Qk], 0);

        //----------------------------------------------------------------------
        // swap rows k and n-1 of L and U
        //----------------------------------------------------------------------
        // swap rows k and n-1 of U           % O(1) time
        v = U->v[k];       U->v[k] = U->v[n-1];    U->v[n-1] = v;
        // column swapping for L has been done by the caller: SPEX_LUU

        // update row permutation to swap rows of L implicitly
        P[k] = Pks;          P[n-1] = Pk;
        P_inv[Pks] = k;   P_inv[Pk] = n-1;          Pks = Pk;
    }
    else // ks < n
    {
        //----------------------------------------------------------------------
        // perform backtracking for frame ks
        //----------------------------------------------------------------------
        // find the scale for backtracking
        if (k == 0)
        {
            // pending_scale = 1/sd(ks-1)
            SPEX_MPQ_SET_UI(pending_scale, 1, 1);
        }
        else
        {
            // pending_scale = sd(k-1)/sd(ks-1)
            SPEX_MPQ_SET_Z(pending_scale, sd[k-1]);
        }
        SPEX_MPQ_SET_DEN(pending_scale, sd[tmp_ks-1]);
        // remove common factor in mpq_den and mpq_num
        SPEX_MPQ_CANONICALIZE(pending_scale);

        // S(:,ks) = pending_scale*S(:,ks)
        SPEX_MPQ_MUL(SL(tmp_ks),
                                SL(tmp_ks), pending_scale);
        SPEX_MPQ_MUL(SU(tmp_ks),
                                SU(tmp_ks), pending_scale);
        // sd(ks) = sd(ks)*pending_scale
        SPEX_MPZ_DIVEXACT(sd[tmp_ks],
                                sd[tmp_ks], SPEX_MPQ_DEN(pending_scale));
        SPEX_MPZ_MUL(sd[tmp_ks],
                                sd[tmp_ks], SPEX_MPQ_NUM(pending_scale));

        //----------------------------------------------------------------------
        // swap rows and columns k and ks of L and U
        //----------------------------------------------------------------------
        // swap entries in sd
        SPEX_MPZ_SWAP(sd[tmp_ks], sd[k]);

        // swap columns k and ks of L        % O(1) time
        v = L->v[k];       L->v[k] = L->v[tmp_ks];    L->v[tmp_ks] = v;
        // swap rows k and ks of U           % O(1) time
        v = U->v[k];       U->v[k] = U->v[tmp_ks];    U->v[tmp_ks] = v;

        // update row/col permutation to swap rows of L and cols of U implicitly
        Q[k] = Qks;          Q[tmp_ks] = Qk;
        Q_inv[Qks] = k;      Q_inv[Qk] = tmp_ks;          Qks = Qk;
        P[k] = Pks;          P[tmp_ks] = Pk;
        P_inv[Pks] = k;      P_inv[Pk] = tmp_ks;          Pks = Pk;

        //----------------------------------------------------------------------
        // update entries in frames between k and ks
        //----------------------------------------------------------------------
        if (tmp_ks > k+1)
        {
            // get the scale for entries between frames k and ks % O(1) time
            // pending_scale = sd(k)/sd (ks);
            SPEX_MPQ_SET_Z(pending_scale, sd[k]);
            SPEX_MPQ_SET_DEN(pending_scale, sd[tmp_ks]);
            SPEX_MPQ_CANONICALIZE(pending_scale);
            // scale entries in frames k+1:ks-1
            for (j = k+1; j < tmp_ks; j++)
            {
                // S(:,j) = S(:,j)*pending_scale;
                SPEX_MPQ_MUL(SL(j), SL(j), pending_scale);
                SPEX_MPQ_MUL(SU(j), SU(j), pending_scale);
                // sd(j) = sd(j)*pending_scale;
                SPEX_MPZ_DIVEXACT(sd[j],
                                        sd[j], SPEX_MPQ_DEN(pending_scale));
                SPEX_MPZ_MUL(sd[j],
                                        sd[j], SPEX_MPQ_NUM(pending_scale));
            }
        }
    }

    //--------------------------------------------------------------------------
    // perform IPGE for row ks, skip IPGE for column since it is all zero
    //--------------------------------------------------------------------------
    int64_t pks, cks, last_nz_b4_ks = k-1;
    // sgn = sgn( S(2,ks) - 1 )
    SPEX_MPQ_CMP_UI(&sgn, SU(tmp_ks), 1, 1);
    // initialize history vector h
    for (pks = 0; pks < Uk_dense_row->nz; pks++)
    {
        cks = Uk_dense_row->i[pks];
        // Lk_dense_col or Uk_dense_row are initialized with no explicit zero
        // for column/row 0 (may contain explicit zero for column/row j>0).  And
        // entries in h are maintained to be >= -1. Therefore, with such
        // initialization, entry with h > -1 is clearly not in nnz pattern and
        // any entry in the nnz pattern with h = -1 must be nonzero.  In all,
        // any explicit zero with h >= -1 must not be in the nnz pattern.
        //
        // REMARK:
        // This is useful only for IPGE update, or when there is possibility of
        // fillin.
        h[cks] = SPEX_FLIP(k-1);

        // no need to apply S(2,ks) if S(2,ks) == 1
        if (sgn == 0) { continue;}

        // apply S(2,tmp_ks) to U(ks,cks) todo: try not to do this
        // This must be done so that the following IPGE update will have the
        // result in integer domain
        SPEX_MPZ_DIVEXACT(Uk_dense_row->x[cks],
                                     Uk_dense_row->x[cks],
                                     SPEX_MPQ_DEN(SU(tmp_ks)));
        SPEX_MPZ_MUL     (Uk_dense_row->x[cks],
                                     Uk_dense_row->x[cks],
                                     SPEX_MPQ_NUM(SU(tmp_ks)));
    }
    // the value of S(2,ks) will be updated at the end of this function

    for (j = k; j < tmp_ks; j++)
    {
        Qj = Q[j];
        // if swapping with the inserted column vk, we handle the first loop
        // j==k seperatedly since the entry has been inserted to L->v[k]
        if (ks == n && j == k)
        {
            // L->v[k] has been updated with the inserted column vk in this
            // case, and there is at most one entry in L->v[k] other than the
            // pivot L(P[k], k) (since this function is called), which is
            // L(P[n-1], k). The pivot entry L(P[k], k) is kept as the first
            // entry in the nnz list, while the second entry (if exists) in the
            // list would be L(P[n-1],k).
            // Therefore, if there is only 1 entry in L->v[k], then
            // L(P[n-1],k) == 0
            if (L->v[k]->nz == 1) { continue; }

            // only need to perform IPGE for U(n-1, Q(n-1)) (i.e., U(k,Q(n-1))
            // before permutation) since there is only one off-diagonal nnz in
            // row k of U (i.e., row n-1 before permutation), which is
            // U(k,Q(n-1)). Since U(k,Q(k)) has not yet been inserted to row k
            // of U, U(k,Q(n-1)) can be found as U->v[k]->x[0], and it has no
            // pending scale, i.e., S(2,k) = 1.
            //
            // U(n-1,Q(n-1)) = (U(n-1, Q(n-1))*L(P[k], k)-
            //                  L(P[n-1], k)*U(k, Q(n-1)))/sd[k-1]
            // NOTE: S(1,k) = 1, so the above result is unscaled
            SPEX_MPZ_MUL(Uk_dense_row->x[Qks],
                                    Uk_dense_row->x[Qks], L->v[k]->x[0]);
            SPEX_MPZ_SUBMUL(Uk_dense_row->x[Qks],
                                    U->v[k]->x[0], L->v[k]->x[1]);
            if (k > 0)
            {
                SPEX_MPZ_DIVEXACT(Uk_dense_row->x[Qks],
                                    Uk_dense_row->x[Qks], sd[k-1]);
            }

            // update history vector
            h[Qk]  = SPEX_UNFLIP(h[Qk]);// U(tmp_ks,Q[k]) is up-to-date
            h[Qks] = SPEX_FLIP(k); // U(tmp_ks,Q[tmp_ks]) is in k-th IPGE

            // update index for last nonzero before ks-th entry
            last_nz_b4_ks = k;

            continue;
        }
        else
        {
            SPEX_MPZ_SGN(&sgn, Uk_dense_row->x[Qj]);
            // skip if U(ks, Q[j]) == 0
            if (sgn == 0) { continue; }
        }

        // perform j-th IPGE update for U(ks,:)
        SPEX_CHECK(spex_update_ipge(Uk_dense_row, h, NULL, U->v[j], Q, Q_inv,
            (const SPEX_matrix)rhos, j));
        // update index for last nonzero before ks-th entry
        last_nz_b4_ks = j;

        // insert new entry L(P(ks), j) to L and swap its value with U(ks, Q(j))
        SPEX_CHECK(spex_update_insert_new_entry(Uk_dense_row->x[Qj], L->v[j],
            SL(j), Pks, option));
        // reset U(ks, Q[j])=0
        SPEX_MPZ_SET_UI(Uk_dense_row->x[Qj], 0);
    }
    if (ks == n)
    {
        SPEX_MPZ_SGN(&sgn, Uk_dense_row->x[Qks]);
        if (sgn == 0)
        {
            // triggered by Tcov/Mats4Tcov/SPEX_Update/mat4.txt, which gives
            // the following frame matrix
            // 1 0 0 1
            // 0 1 0 1
            // 0 0 1 0
            // 0 0 0 1
            // ^
            // |
            // update this column with [1; 0; 0; 1]
            // run the following in Tcov folder for more details:
            // ./tcov_test 0 1 Mats4Tcov/SPEX_Update/mat4.txt
            SPEX_FREE_ALL;
            return SPEX_SINGULAR;
        }
    }

    // There must be at least one nonzero in U(ks, Q(k:ks-1)) for this case.
    // Otherwise, this will be handled by spex_dppu1(...). Therefore,
    // last_nz_b4_ks must have updated at least once and be greater than k-1
    ASSERT(last_nz_b4_ks > k-1);

    if (ks == n)
    {
        // There should be only 1 entry in row tmp_ks (i.e., n-1) of U, which
        // is U(ks,Q[tmp_ks]). History update will be performed right after (by
        // checking last_nz_b4_ks). In addition, the history vector for all
        // other entries in row tmp_ks are up-to-date (i.e., maintained as
        // >=-1)
        h[Qks] = SPEX_UNFLIP(h[Qks]);
        last_nz_b4_ks = h[Qks];
    }
    else
    {
        // perform history update up to (last_nz_b4_ks)-th IPGE iteration
        // and remove zero from row ks of U
        *jnext = n;
        for (pks = 0; pks < Uk_dense_row->nz;)
        {
            // column index in row ks of U
            cks = Uk_dense_row->i[pks];
            h[cks] = SPEX_UNFLIP(h[cks]);
            int64_t real_cks = Q_inv[cks];
            SPEX_MPZ_SGN(&sgn, Uk_dense_row->x[cks]);
            if (sgn == 0 || real_cks < tmp_ks)
            {
                // Remove indices of explicit zeros from nnz pattern. In
                // addition, all entries in U(ks, Q(k:ks-1)) should be zero
                /*if (sgn != 0)
                {
                    SPEX_MPZ_SET_UI(Uk_dense_row->x[cks], 0);
                }*/
                Uk_dense_row->nz--;
                Uk_dense_row->i[pks] = Uk_dense_row->i[Uk_dense_row->nz];
                continue;
            }


            // update the index of next off-diagonal nnz entry
            if (real_cks > tmp_ks && real_cks < *jnext)
            {
                *jnext = real_cks;
            }

            if (h[cks] < last_nz_b4_ks) // require history update
            {
                // U(ks,cks) = (U(ks,cks)*sd(last_nz_b4_ks))/sd(h[cks]);
                SPEX_MPZ_MUL(Uk_dense_row->x[cks],
                                      Uk_dense_row->x[cks], sd[last_nz_b4_ks]);
                if (h[cks] > -1)
                {
                    SPEX_MPZ_DIVEXACT(Uk_dense_row->x[cks],
                                     Uk_dense_row->x[cks], sd[h[cks]]);
                }
            }
            pks++;
        }
    }

    if (last_nz_b4_ks != tmp_ks-1)
    {
        // update S(2,ks) = sd(ks-1)/sd(last_nz_b4_ks)
        SPEX_MPQ_SET_Z(SU(tmp_ks), sd[tmp_ks-1]);
        SPEX_MPQ_SET_DEN(SU(tmp_ks), sd[last_nz_b4_ks]);
        SPEX_MPQ_CANONICALIZE(SU(tmp_ks));
        // sd(ks) = U(ks,Q(ks))*S(3,ks)
        SPEX_MPZ_DIVEXACT(sd[tmp_ks],Uk_dense_row->x[Qks],
                                     SPEX_MPQ_DEN(SU(tmp_ks)));
        SPEX_MPZ_MUL(sd[tmp_ks], sd[tmp_ks],
                                     SPEX_MPQ_NUM(SU(tmp_ks)));
    }
    else
    {
        // update S(2,ks) = 1
        SPEX_MPQ_SET_UI(SU(tmp_ks), 1, 1);
        // sd(ks) = U(ks,Q[ks])
        SPEX_MPZ_SET(sd[tmp_ks], Uk_dense_row->x[Qks]);
    }


    if (ks == n)
    {
        // move data from Uk_dense_row, there is only one entry that needs
        // to move, which is U(k,Q[n-1])
        // set U(n-1,n-1)=L(n-1,n-1)=Uk_dense_row[Q[n-1]]
        SPEX_MPZ_SWAP(U->v[n-1]->x[0], Uk_dense_row->x[Qks]);
        SPEX_MPZ_SET (L->v[n-1]->x[0], U->v[n-1]->x[0]     );
        U->v[n-1]->i[0] = Qks;
        U->v[n-1]->nz = 1;
        L->v[n-1]->i[0] = Pks;
        // S(1,n-1) = S(2, n-1)
        SPEX_MPQ_SET(SL(tmp_ks), SU(n-1));
    }
#ifdef SPEX_DEBUG
    // no need to update L(P(ks), ks), which will be updated in the last loop
    // L(P[ks],ks) = U(ks,Q[ks]);
    // only do so in the debug mode for simpler verification
    else
    {
        SPEX_MPZ_SET(Lk_dense_col->x[Pks],
                                Uk_dense_row->x[Qks]);
        SPEX_MPQ_SET(SL(tmp_ks), SU(tmp_ks));
    }
#endif

    SPEX_FREE_ALL;
    return SPEX_OK;
}
