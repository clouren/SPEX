//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_dppu1.c: perform diagonal permutation pivot update
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

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
    SPEX_MPZ_CLEAR(Lksk);            \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPZ_CLEAR(tmpz);

#include "spex_lu_update_internal.h"

SPEX_info spex_dppu1
(
    SPEX_mat *L,  // matrix L
    SPEX_mat *U,  // matrix U
    mpq_t *S,        // array of size 3*n that stores pending scales
    mpz_t *d,        // array of size n that stores the unscaled pivot
    mpz_t *sd,       // array of size n that stores the scaled pivot
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    int64_t *inext,  // the index of first off-diag entry in col k of L
    int64_t *h,      // allocated vector that can be used for history vector.
                     // All entries are maintained to be >= -1
    int64_t *Q,      // column permutation
    int64_t *Q_inv,  // inverse of column permutation
    int64_t *P,      // row permutation
    int64_t *P_inv,  // inverse of row permutation
    int64_t *Ldiag,  // L(k,k) can be found as L->v[k]->x[Ldiag[k]]
    const int64_t *Uci,// the row index for col-wise nnz pattern of U
    const int64_t *Ucp,// col pointers for col-wise nnz pattern of U
    const int64_t *Ucx,// the value of k-th entry is found as
                       // U->v[Uci[k]]->x[Ucx[k]]
    const int64_t k,   // current column index 0 <= k < n
    const int64_t ks   // index of the diagonal to be swapped with, [0,n)
)
{
    printf("using dppu1 swapping k(%ld) and ks(%ld)\n",k,ks);
    // initialize workspace
    SPEX_info info;
    int sgn, Lksk_sgn;
    int64_t pk, ck, pks, cks, tmpi, j, n = U->n;
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
    // NOTE: U(k,Q(ks)) must be nnz, otherwise it would cause singularity, since
    //       U(k,Q(k+1:n)) will be all zeros
    // -------------------------------------------------------------------------
    if (ks == n)
    {
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
        SPEX_CHECK(SPEX_mpz_sgn(&(Lksk_sgn), Lk_dense_col->x[P[n-1]]));
        if (Lksk_sgn == 0)
        {
            SPEX_CHECK(SPEX_mpz_divexact(U->v[n-1]->x[0],
                                         U->v[n-1]->x[0], sd[n-2]));
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_fdiv_q(U->v[n-1]->x[0],
                                       U->v[n-1]->x[0], sd[n-2]));
            // pending_scale = S(1,k)*S(3,k)
            SPEX_CHECK(SPEX_mpq_mul(pending_scale, SPEX_LUU_2D(S, 1, k),
                                    SPEX_LUU_2D(S, 3, k) ));
            // L(P(n-1),k) = L(P(n-1),k)*pending_scale, which should be integer
            SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[P[n-1]],
                       Lk_dense_col->x[P[n-1]], SPEX_MPQ_DEN(pending_scale)));
            SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[P[n-1]],
                       Lk_dense_col->x[P[n-1]], SPEX_MPQ_NUM(pending_scale)));
            // tmpz = ceil(U(k,Q(n-1))*L(P(n-1),k)/U(k,Q(k))
            SPEX_CHECK(SPEX_mpz_mul(tmpz, Uk_dense_row->x[Q[n-1]],
                                    Lk_dense_col->x[P[n-1]]));
            SPEX_CHECK(SPEX_mpz_cdiv_q(tmpz, tmpz, Uk_dense_row->x[Q[k]]));
            // U(n-1,Q(n-1)) = U(n-1,Q(n-1))+tmpz
            SPEX_CHECK(SPEX_mpz_add(U->v[n-1]->x[0], U->v[n-1]->x[0], tmpz));
        }

        // ---------------------------------------------------------------------
        // scale entries in frames k+1:n-2
        // ---------------------------------------------------------------------
        // since the value in Uk_dense_row[Q[k]] will not be used, we use it to
        // hold the original value of sd[k] before swapping columns and rows of
        // k and n-1. Then we set sd[k] to d[k]
        SPEX_CHECK(SPEX_mpz_swap(Uk_dense_row->x[Q[k]], sd[k]));
        SPEX_CHECK(SPEX_mpz_set(sd[k], d[k]));

        if (n > k+2) // n-1 > k+1
        {
            // pending_scale = sd(k)/Uk_dense_row[Q[k]]
            SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k]));
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, Uk_dense_row->x[Q[k]]));
            SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

            for (j = k+1; j < n-1; j++)
            {
                // S(3,k+1:n-2) = S(3,k+1:n-2)*pending_scale;
                SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 3, j),
                                        SPEX_LUU_2D(S, 3, j), pending_scale));
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
        tmpi = P[k];       P[k] = P[n-1];          P[n-1] = tmpi;
        P_inv[P[k]] = k;   P_inv[tmpi] = n-1;

        // move data from Uk_dense_row, there is only one entry that needs
        // to move, which is U(k,Q[n-1])
        // set U(n-1,Q(n-1))=L(P(n-1),n-1)=Uk_dense_row[Q[n-1]]
        SPEX_CHECK(SPEX_mpz_swap(U->v[n-1]->x[0], Uk_dense_row->x[Q[n-1]]));
        SPEX_CHECK(SPEX_mpz_set (L->v[n-1]->x[0], U->v[n-1]->x[0]     ));
        U->v[n-1]->i[0] = Q[n-1];
        U->v[n-1]->nz = 1;
        L->v[n-1]->i[0] = P[n-1];
        L->v[n-1]->nz = 1;

        // ---------------------------------------------------------------------
        // set d, sd and S for frame n-1
        // ---------------------------------------------------------------------
        // the scaling factor for frame n-1 same as that for row k of U, while
        // multiply with the scaling factor due to the IPGE update.
        SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 3, n-1),
                                SPEX_LUU_2D(S, 3, k), SPEX_LUU_2D(S, 2, k)));
        // get the scale for IPGE update: pending_scale = sd(n-2)/sd(k-1);
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[n-2]));
        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k-1]));
            SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        }
        SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 3, n-1),
                                SPEX_LUU_2D(S, 3, n-1), pending_scale));
        SPEX_CHECK(SPEX_mpq_set_ui(SPEX_LUU_2D(S, 1, n-1), 1, 1));
        SPEX_CHECK(SPEX_mpq_set_ui(SPEX_LUU_2D(S, 2, n-1), 1, 1));

        // d[n-1] = U(n-1,Q(n-1))
        SPEX_CHECK(SPEX_mpz_set(d[n-1],  U->v[n-1]->x[0]));
        SPEX_CHECK(SPEX_mpz_divexact(sd[n-1],
                                 d[n-1], SPEX_MPQ_DEN(SPEX_LUU_2D(S, 3, n-1))));
        SPEX_CHECK(SPEX_mpz_mul(sd[n-1],
                                sd[n-1], SPEX_MPQ_NUM(SPEX_LUU_2D(S, 3, n-1))));

        SPEX_FREE_ALL;
        return SPEX_OK;
    }

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

    // S(3, ks) = pending_scale*S(3,ks)
    SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 3, ks),
                            SPEX_LUU_2D(S, 3, ks), pending_scale));

    SPEX_CHECK(SPEX_mpz_sgn(&(Lksk_sgn), Lk_dense_col->x[P[ks]]));
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
        // update the scale for col ks of L due to backtracking
        // S(1,ks) = S(1,ks)*S(3,ks) = S(1,ks)*S(3,ks)*sd(k-1)/sd(ks-1)
        SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 1, ks),
                                SPEX_LUU_2D(S, 1, ks), SPEX_LUU_2D(S, 3, ks)));

        // S(2,ks) = S(2,ks)*S(3,ks) = S(2,ks)*S(3,ks)*sd(k-1)/sd(ks-1)
        SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 2, ks),
                                SPEX_LUU_2D(S, 2, ks), SPEX_LUU_2D(S, 3, ks)));
        // S(3, ks) = 1
        SPEX_CHECK(SPEX_mpq_set_ui(SPEX_LUU_2D(S, 3, ks), 1, 1));

        // pending_scale = S(1,k)*S(3,k)
        SPEX_CHECK(SPEX_mpq_mul(pending_scale,
                                SPEX_LUU_2D(S, 1, k), SPEX_LUU_2D(S, 3, k) ));
        // assign Lksk = L(P(ks),k)*pending_scale, which should be integer
        SPEX_CHECK(SPEX_mpz_divexact(Lksk, Lk_dense_col->x[P[ks]],
                                     SPEX_MPQ_DEN(pending_scale)));
        SPEX_CHECK(SPEX_mpz_mul(Lksk, Lksk, SPEX_MPQ_NUM(pending_scale)));

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
        h[Q[k]] = -1;
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
                SPEX_CHECK(SPEX_mpz_cdiv_q(tmpz, tmpz, Uk_dense_row->x[Q[k]]));
                // U(ks,cks) = floor(U(ks,cks)*S(2,ks))
                SPEX_CHECK(SPEX_mpz_mul(U->v[ks]->x[pks], U->v[ks]->x[pks],
                                        SPEX_MPQ_NUM(SPEX_LUU_2D(S, 2, ks))));
                SPEX_CHECK(SPEX_mpz_fdiv_q(U->v[ks]->x[pks], U->v[ks]->x[pks],
                                        SPEX_MPQ_DEN(SPEX_LUU_2D(S, 2, ks))));
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
                                        SPEX_MPQ_DEN(SPEX_LUU_2D(S, 2, ks))));
                SPEX_CHECK(SPEX_mpz_mul(U->v[ks]->x[pks], U->v[ks]->x[pks],
                                        SPEX_MPQ_NUM(SPEX_LUU_2D(S, 2, ks))));

                // update sd[ks] = U(ks, Q(ks))
                if (cks == Q[ks])
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
                SPEX_CHECK(SPEX_vector_realloc(U->v[ks], U->v[ks]->nz+count));
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
                        U->v[ks]->x[pks], Uk_dense_row->x[Q[k]]));
                    pks++;

                    h[ck] = -1;
                    count--;
                }
            }
            U->v[ks]->nz = pks;
        }

        // d[ks] = L(P(ks), ks)
        SPEX_CHECK(SPEX_mpz_set(d[ks], L->v[ks]->x[Ldiag[ks]]));
        // reset S(2,ks)
        SPEX_CHECK(SPEX_mpq_set_ui(SPEX_LUU_2D(S, 2, ks), 1, 1));

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
    tmpi = Q[k];       Q[k] = Q[ks];          Q[ks] = tmpi;
    Q_inv[Q[k]] = k;   Q_inv[tmpi] = ks;
    tmpi = P[k];       P[k] = P[ks];          P[ks] = tmpi;
    P_inv[P[k]] = k;   P_inv[tmpi] = ks;

    // update Ldiag[k] = Ldiag[ks]
    Ldiag[k] = Ldiag[ks];

    // swap entries in d, sd and S
    SPEX_CHECK(SPEX_mpz_swap(sd[k], sd[ks]));
    SPEX_CHECK(SPEX_mpz_swap( d[k],  d[ks]));
    SPEX_CHECK(SPEX_mpq_swap(SPEX_LUU_2D(S, 1, ks), SPEX_LUU_2D(S, 1, k)));
    SPEX_CHECK(SPEX_mpq_swap(SPEX_LUU_2D(S, 2, ks), SPEX_LUU_2D(S, 2, k)));
    SPEX_CHECK(SPEX_mpq_swap(SPEX_LUU_2D(S, 3, ks), SPEX_LUU_2D(S, 3, k)));

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
            // S(3,k+1:ks-1) = S(3,k+1:ks-1)*pending_scale;
            SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 3, j),
                                    SPEX_LUU_2D(S, 3, j), pending_scale));
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
    if (Lksk_sgn == 0 || ks == n-1)
    {
        // get the scale for IPGE update: pending_scale = sd(ks-1)/sd (k-1);
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[ks-1]));
        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k-1]));
            SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        }

        // S(3, ks) = S(3, ks)*pending_scale
        SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 3, ks),
                                SPEX_LUU_2D(S, 3, ks), pending_scale));
        // sd(ks) = sd(ks)*pending_scale;
        SPEX_CHECK(SPEX_mpz_divexact(sd[ks], sd[ks],
                                     SPEX_MPQ_DEN(pending_scale)));
        SPEX_CHECK(SPEX_mpz_mul(sd[ks], sd[ks], SPEX_MPQ_NUM(pending_scale)));
    }
    else
    {
        // skip scaling for U for 1 IPGE iteration
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k]));
        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k-1]));
            SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        }
        // S(2,ks) = S(2,ks)*pending_scale
        SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 2, ks),
                                SPEX_LUU_2D(S, 2, ks), pending_scale));

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
        ASSERT(pks == Lk_dense_col->nz);

        // NOTE: This will cause fillin in the ks(th) column of L,
        //       Since there is no subset relation between nnz pattern in
        //       L(:,ks) and L(:, k). Both could have explicit zero(s).
        //       L(:,k) can be jumbled.
        *inext = n;
        for (pk = 0; pk < L->v[k]->nz; pk++)
        {
            // row index in column k of L
            ck = L->v[k]->i[pk];
            // exclude L(P[k],k)
            if (ck == P[k])
            {
                continue;
            }

            // L(ck,ks) = (L(ck,ks)*d(k)-L(ck,k)*L(P(k),ks))/sd(k-1);
            // use L(P(k), ks) without applying any pending scale factor to
            // keep it in the same pending scale factor as the rest of entries
            // in column k of L, so that the skipped scaling for col k of L can
            // still be skipped when performing IPGE update
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[ck]));
            if (sgn != 0) // L(ck, ks) != 0
            {
                SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ck],
                                        Lk_dense_col->x[ck], d[k]));
            }
            else if (h[ck] >= -1) // this entry was not in nnz pattern
            {
                // insert new entry in the nonzero pattern
                ASSERT(pks <= Lk_dense_col->nzmax);
                Lk_dense_col->i[pks] = ck;
                pks++;
            }
            SPEX_CHECK(SPEX_mpz_submul(Lk_dense_col->x[ck],
                                       Lk_dense_col->x[P[k]], L->v[k]->x[pk]));
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[ck]));
            if (sgn != 0) // L(ck, ks) != 0
            {
                if (k > 0)
                {
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ck],
                                                 Lk_dense_col->x[ck], sd[k-1]));
                }

                // check if this will be the first off diagonal entry in L(P,ks)
                if (P_inv[ck] < *inext)
                {
                    ASSERT(P_inv[ck] > ks);
                    // inext is the row index of the found first off-diagonal
                    // entry in L(P,ks)
                    *inext = P_inv[ck];
                }
            }
            
            // update h[ck] to mark Lk_dense_col[ck] need no further update
            h[ck] = -1;
        }
        // update the number of nnz with pks-1 since L(P(k),ks) will be deleted
        Lk_dense_col->nz = pks-1;

        for (pks = 0; pks < Lk_dense_col->nz; pks ++)
        {
            // row index in column ks of L
            cks = Lk_dense_col->i[pks];

            if (h[cks] < -1) //only need to update entries that were not updated
            {
                // reset history vector to any value >= -1
                h[cks] = -1;
                // L(P(k), ks) should be removed from nnz pattern
                if (cks == P[k])
                {
                    // move the row index of last nonzero to current position
                    Lk_dense_col->i[pks] = Lk_dense_col->i[Lk_dense_col->nz];
                    SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col->x[cks], 0));
                    continue;
                }

                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[cks]));
                if (sgn != 0) // L(cks, ks) != 0
                {
                    // L(cks,ks) = (L(cks,ks)*d(k))/sd(k-1);
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[cks],
                                            Lk_dense_col->x[cks], d[k]));
                    if (k > 0)
                    {
                        SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[cks],
                                            Lk_dense_col->x[cks], sd[k-1]));
                    }
                    // check if this will be the 1st off-diag entry in L(P,ks)
                    if (P_inv[cks] < *inext && P_inv[cks] > ks)
                    {
                        // inext is the row index of the found first
                        // off-diagonal entry in L(P,ks)
                        *inext = P_inv[cks];
                    }
                }
            }
        }

        // S(1,ks) = S(1,ks)*S(1,k)
        SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 1, ks),
                                SPEX_LUU_2D(S, 1, ks), SPEX_LUU_2D(S, 1, k)));
        
        // skip the rest of IPGE iterations
        // pending_scale = sd(ks-1)/sd(k);
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[ks-1]));
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[k]));
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));
        // S(3, ks) = S(3, ks)*pending_scale;
        SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 3, ks),
                                SPEX_LUU_2D(S, 3, ks), pending_scale));
        // sd(ks) = L(P(ks),ks)*S(1,ks)*S(3,ks);
        SPEX_CHECK(SPEX_mpq_mul(pending_scale,
                                SPEX_LUU_2D(S, 3, ks), SPEX_LUU_2D(S, 1, ks)));
        SPEX_CHECK(SPEX_mpz_divexact(sd[ks], Lk_dense_col->x[P[ks]],
                                SPEX_MPQ_DEN(pending_scale)));
        SPEX_CHECK(SPEX_mpz_mul(sd[ks], sd[ks],
                                SPEX_MPQ_NUM(pending_scale)));
    }
    SPEX_FREE_ALL;
    return SPEX_OK;
}
