//------------------------------------------------------------------------------
// SPEX_Update/spex_update_cppu.c: perform column permutation pivot update
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is called to perform column permutation pivot update
// when the submatrix (formed by rows and columns k to ks) has the following
// pattern
//       x 0 0 0 x       <- row k
//       . x . . .
//       . . x . .
//       . . . x .
//       . . . . x       <- row ks
//
//       ^       ^
//       |       |
//     col k   col ks
//       
// This function will swap columns k and ks in L and U. Noted that the columns
// of U are permuted implicitly via the permutation matrix based on Q.

#define SPEX_FREE_ALL                \
    SPEX_MPZ_CLEAR(Uiks);            \
    SPEX_MPQ_CLEAR(one);             \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPZ_CLEAR(tmpz);

#include "spex_update_internal.h"

#define SL(k) (L->v[(k)]->scale)
#define SU(k) (U->v[(k)]->scale)

SPEX_info spex_update_cppu
(
    SPEX_matrix *L,     // matrix L
    SPEX_matrix *U,     // matrix U
    SPEX_matrix *rhos,// array of scaled pivots
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    int64_t *inext,  // the index of first off-diag entry in col k of L
    int64_t *jnext,  // the index of first off-diag entry in row k of U
    int64_t *h,      // allocated vector that can be used for history vector.
                     // All entries are maintained to be >= -1
    int64_t *Q,      // column permutation
    int64_t *Q_inv,  // inverse of column permutation
    int64_t *P,      // row permutation (unchanged on output)
    int64_t *P_inv,  // inverse of row permutation (unchanged on output)

    // the col-wise nnz pattern of U can be NULL when ks == n
    const int64_t *Uci_ks,// the row index for nnz pattern of U(k+1:ks-1,Q[ks])
    const int64_t *Ucx_ks,// the value of i-th entry is found as
                     // U->v[Uci[i]]->x[Ucx[i]]
    const int64_t Uc_ks_nz,// # of nnz in U(k+1:ks-1,Q[ks])
    const int64_t k, // current column index 0 <= k < n
    const int64_t ks // index of the diagonal to be swapped with, [0,n)
)
{
#ifdef SPEX_DEBUG
    printf("using cppu swapping k(%ld) and ks(%ld), Uc_nz=%ld\n",k,ks,Uc_ks_nz);
#endif
    // initialize workspace
    SPEX_info info;
    int sgn, r;
    int64_t pk, ck, pks, cks, pi, ci, i, j, n = U->n;
    int64_t Lk_nz = 0;
    // the pointer for U(k,Q(ks)) = U->v[k]->x[Ucx[Ucp_k_ks]]
    int64_t Ucp_k_ks = -1;
    int64_t Qk = Q[k], Qks, Pk = P[k], Pks, Pi;
    *inext = n;
    *jnext = n;
    mpz_t *sd = rhos->x.mpz;

    mpq_t pending_scale, one;
    SPEX_MPQ_SET_NULL(pending_scale);
    mpz_t Uiks, tmpz; SPEX_MPZ_SET_NULL(Uiks); SPEX_MPZ_SET_NULL(tmpz);
    SPEX_CHECK(SPEX_mpq_init(pending_scale));
    SPEX_CHECK(SPEX_mpq_init(one));
    SPEX_CHECK(SPEX_mpq_set_ui(one, 1, 1));
    SPEX_CHECK(SPEX_mpz_init(Uiks));
    SPEX_CHECK(SPEX_mpz_init(tmpz));

    if (ks == n)
    {
#ifdef SPEX_DEBUG
        printf("ks==n using cppu\n");
#endif
        Qks = Q[n-1];
        // since the value in Uk_dense_row[Q[k]] will not be used, we use it to
        // hold the original value of sd[k] before swapping column k with
        // column n-1. Then we set sd[k] to the new pivot of column k of L
        // vk[P[k]], which is kept as first entry in the nnz list.
        SPEX_CHECK(SPEX_mpz_swap(Uk_dense_row->x[Qk], sd[k]));
        SPEX_CHECK(SPEX_mpz_set(sd[k], L->v[k]->x[0]));

        // get the scale for entries between frames k and n-1
        // pending_scale = sd(k)/Uk_dense_row[Q[k]]
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k]));
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, Uk_dense_row->x[Qk]));
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

        // if the inserted column is used, we won't need to perform
        // backtracking. Instead, we just need to perform RwSOP.
        //
        // If U(k, Q[n-1]) == 0, RwSOP is simplified as pure scaling for all
        // frame from k+1:n-1. Otherwise (which won't happen due to heuristic),
        // we need to first compute the (n-1)-th IPGE iteration for the column
        // k, and use the result to perform RwSOP on column n-1 of U.
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row->x[Qks]));
        if (sgn == 0)
        {
            // just need to perform RwSOP by scaling all frame k+1:n-1
            for (j = k+1; j < n; j++)
            {
                // S(:,k+1:n-1) = S(:,k+1:n-1)*pending_scale;
                SPEX_CHECK(SPEX_mpq_mul(SL(j), SL(j), pending_scale));
                SPEX_CHECK(SPEX_mpq_mul(SU(j), SU(j), pending_scale));
                // sd(k+1:n-1) = sd(k+1:n-1)*pending_scale;
                SPEX_CHECK(SPEX_mpz_divexact(sd[j],
                                        sd[j], SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(sd[j],
                                        sd[j], SPEX_MPQ_NUM(pending_scale)));
            }

            // move data from Uk_dense_row, but there is no entry that needs
            // to move
            U->v[k]->nz = 0;
        }
        else
        {
            // apply pending scale to row k of U, which has 1 nnz U(k,Q[n-1])
            // U(k,Q[n-1]) = U(k, Q[n-1])*S(2,k)
            SPEX_CHECK(SPEX_mpz_divexact(Uk_dense_row->x[Qks],
                      Uk_dense_row->x[Qks], SPEX_MPQ_DEN(SU(k))));
            SPEX_CHECK(SPEX_mpz_mul(Uk_dense_row->x[Qks],
                      Uk_dense_row->x[Qks], SPEX_MPQ_NUM(SU(k))));
            // S(:,k) will be set to 1 after calling this function in SPEX_LUU

            // perform 1 IPGE iteration on Lk_dense_col using vk (which has
            // been inserted as L->v[k]) and update the history vector. Then
            // use L to perform the remaining IPGE update till (n-1)-th
            // iteration. Finally, use the result to update column n-1 of U.
            //
            // It should be noted that, since the resulted column in the
            // (n-1)-th IPGE iteration is computed using column k of L instead
            // of the inserted column, its sign should be flipped when applying
            // RwSOP.
            if (Lk_dense_col->nz == 1)
            {
                SPEX_CHECK(SPEX_mpz_set(Lk_dense_col->x[Pk], sd[k]));
                h[Pk] = SPEX_FLIP(k-1);
            }
            else
            {
                SPEX_CHECK(SPEX_mpz_gcd(tmpz, sd[k-1], SPEX_MPQ_NUM(SL(k))));
                for (pk = 0; pk < Lk_dense_col->nz; pk++)
                {
                    ck = Lk_dense_col->i[pk];
                    // initialize history vector
                    h[ck] = SPEX_FLIP(k-1);
                    // Instead of apply the real pending scale to Lk_dense_col
                    // before using IPGE in order to maintain the result of
                    // IPGE in the integer domain, we just multiply each entry
                    // with tmpz, which will do the same work. In addition,
                    // there is no need to update the pending scale for
                    // Lk_dense_col, since it won't be needed nor used beyond
                    // this point and Lk_dense_col will be cleared at the end
                    // of LU update process.
                    //SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ck],
                    //          Lk_dense_col->x[ck], SPEX_MPQ_DEN(SL(k))));
                    //SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ck],
                    //          Lk_dense_col->x[ck], SPEX_MPQ_NUM(SL(k))));
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ck],
                                            Lk_dense_col->x[ck], tmpz));
                }
            }

            // skip updating L(P[k],k), so start from pk = 1
            for (pk = 1; pk < L->v[k]->nz; pk++)
            {
                ck = L->v[k]->i[pk];
                
#ifdef SPEX_DEBUG
                // all explicit zeros are removed when put into L->v[k]
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[k]->x[pk]));
                ASSERT (sgn != 0);
#endif

                // L(ck,k) = (L(ck, k)*vk(P[k])-L(P[k],k)*vk(ck))/sd[k-1]
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[ck]));
                if (sgn != 0)
                {
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ck],
                                            Lk_dense_col->x[ck], sd[k]));
                }
                else if (h[ck] >= -1) //this entry wasn't in the nnz pattern
                {
                    // insert new entry in the nonzero pattern
                    Lk_dense_col->i[Lk_dense_col->nz] = ck;
                    Lk_dense_col->nz++;
                }
                SPEX_CHECK(SPEX_mpz_submul(Lk_dense_col->x[ck],
                                           Lk_dense_col->x[Pk],
                                           L->v[k]->x[pk]));
                if (k > 0)
                {
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ck],
                                                 Lk_dense_col->x[ck], sd[k-1]));
                }
                h[ck] = SPEX_FLIP(k);
            }

            // perform IPGE and RwSOP
            for (i = k+1; i < n; i++)
            {
                Pi = P[i];
                // S(:,i) = S(:, i)*pending_scale;
                SPEX_CHECK(SPEX_mpq_mul(SL(i), SL(i), pending_scale));
                SPEX_CHECK(SPEX_mpq_mul(SU(i), SU(i), pending_scale));

                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[Pi]));
                // we compute sd[n-1] differently
                if (sgn == 0 || i != n-1)
                {
                    // sd[i] = sd[i]*pending_scale
                    SPEX_CHECK(SPEX_mpz_mul(sd[i], sd[i],
                                            SPEX_MPQ_NUM(pending_scale)));
                    SPEX_CHECK(SPEX_mpz_divexact(sd[i], sd[i],
                                            SPEX_MPQ_DEN(pending_scale)));

                    // skip if Lk_dense_col[P[i]] == 0
                    if (sgn == 0)       {continue;}
                }

                // perform i-th IPGE update for Lk_dense_col
                SPEX_CHECK(spex_update_ipge(Lk_dense_col, h, NULL, L->v[i],
                    P, P_inv, (const SPEX_matrix*)rhos, i));

                // perform RwSOP for row i with flipped-sign entries in
                // Lk_dense_col. All entries in row i of U must be SCALEUP such
                // that S(:,i)=[S(1,i);1]
                int64_t p = 0; // the pointer to U(i,Q(n-1))
                if ( i != n-1 )
                {
                    p = -1;
                    // set U(i, Q(i)) = sd[i]
                    SPEX_CHECK(SPEX_mpz_set(U->v[i]->x[0], sd[i]));

                    // iterate all nnz in row i of U
                    for (pi = 1; pi < U->v[i]->nz; pi++)
                    {
                        ci = U->v[i]->i[pi];
                        if (ci == Qks)
                        {
                            p = pi;
                        }
                        else //if (Q_inv[ci] < n-1)
                        {
                            // apply S(2,i) to U(i,Q(i:n-2))
                            SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[pi],
                                U->v[i]->x[pi], SPEX_MPQ_DEN(SU(i))));
                            SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[pi],
                                U->v[i]->x[pi], SPEX_MPQ_NUM(SU(i))));
                        }
                    }
                }
                // perform RwSOP to U(i,Q(n-1)), POSSIBLE FILLIN
                // sign is changed here due to column swap
                // U(i,Qks)= U(i,Qks)*S(2,i) +
                //        Lk_dense_col(P[i])*U(k,Qks)/Lk_dense_col[P[k]]
                if (p > -1)
                {
                    SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[p],
                               U->v[i]->x[p], SPEX_MPQ_NUM(SU(i))));
                    SPEX_CHECK(SPEX_mpz_fdiv_q(U->v[i]->x[p],
                               U->v[i]->x[p], SPEX_MPQ_DEN(SU(i))));
                    SPEX_CHECK(SPEX_mpz_mul(tmpz,
                               Lk_dense_col->x[Pi], Uk_dense_row->x[Qks]));
                    SPEX_CHECK(SPEX_mpz_cdiv_q(tmpz, tmpz,
                               Lk_dense_col->x[Pk]));
                    SPEX_CHECK(SPEX_mpz_add(U->v[i]->x[p],U->v[i]->x[p], tmpz));
                }
                else // U(i,Q(n-1)) was not in the nnz pattern
                {
                    p = U->v[i]->nz;
                    // reallocate the nonzero pattern if needed
                    if (p == U->v[i]->nzmax)
                    {
                        SPEX_CHECK(SPEX_vector_realloc(U->v[i],
                            SPEX_MIN(n, 2*(U->v[i]->nzmax)), NULL));
                    }
                    // insert new entry in the nonzero pattern
                    U->v[i]->i[p] = Qks;
                    U->v[i]->nz++;

                    SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[p],
                                  Lk_dense_col->x[Pi], Uk_dense_row->x[Qks]));
                    SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[p],
                                  U->v[i]->x[p], Lk_dense_col->x[Pk]));
                }
                SPEX_CHECK(SPEX_mpq_set_ui(SU(i), 1, 1));
                
                if (i == n-1)
                {
                    // make sure U(n-1,Q(n-1)) != 0
                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[n-1]->x[0]));
                    if (sgn == 0)
                    {
                        SPEX_FREE_ALL;
                        return SPEX_SINGULAR;
                    }

                    // L(P(n-1),n-1) = sd(n-1) = U(n-1, Q(n-1))
                    SPEX_CHECK(SPEX_mpz_set(L->v[n-1]->x[0], U->v[n-1]->x[0]));
                    SPEX_CHECK(SPEX_mpz_set(  sd[n-1]      , U->v[n-1]->x[0]));
                    // S(:,n-1) = ones
                    SPEX_CHECK(SPEX_mpq_set(SL(n-1), one));
                    SPEX_CHECK(SPEX_mpq_set(SU(n-1), one));
                }
            }

            // move data from Uk_dense_row, there is only one entry that needs
            // to move, which is U(k,Q[n-1]) and has no pending scale
            SPEX_CHECK(SPEX_mpz_swap(U->v[k]->x[0], Uk_dense_row->x[Qks]));
            U->v[k]->i[0] = Qks;
            U->v[k]->nz = 1;
        }

        // reset Uk_dense_row->x[Q[k]]=0
        SPEX_CHECK(SPEX_mpz_set_ui(Uk_dense_row->x[Qk], 0));
        SPEX_FREE_ALL;
        return SPEX_OK;
    }
    //-------------------------------------------------------------------------
    // Backtracking column ks of L and U, the backtracking result will be moved
    // to column k of L.
    //-------------------------------------------------------------------------
    // When U(k+1:ks-1,Q(ks)) are all zero(s), there is no need to make a copy
    // of Lk_dense_col. Instead, we can direct operate on Lk_dense_col.
    // Otherwise, Lk_dense_col will firstly be moved to L->v[k]->x in a
    // compressed-column form. Then Lk_dense_col will make a copy of ks-th
    // column of scaled L and scaled U(k:ks-1,Q(ks)). Backtracking will be
    // performed on Lk_dense_col for each nonzero in U(k:ks-1, Q(ks)).
    //-------------------------------------------------------------------------
    Qks = Q[ks];
    Pks = P[ks];

    if (Uc_ks_nz == 0)   // If U(k+1:ks-1,Q(ks)) are all zero(s).
    {
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . 
        // backtracking jumbled sparse column ks of L using 'dense' column k of
        // L and store the result in Lk_dense_col. This will introduce new
        // entry to Lk_dense_col.
        //
        // When explicit zeros in L(:, ks) resulted from exact
        // cancellation in IPGE update were not removed (SLIP LU keeps those
        // zeros in output L and U), nonzero pattern of L(P(ks:n+1),k) should
        // be a subset of L(:,ks). Therefore, the backtracking will need to 
        // simply iterate all nonzero in the L(:,ks), and the final Lk_dense_col
        // will have mostly the same nnz pattern as L(:,ks), except L(P(k),k).
        //
        // However, since we assume the IPGE update results in exact
        // cancellation and the resulted zero is removed from L or U, the
        // subset relation of the nonzero pattern no longer holds. In this
        // case, we need to first update Lk_dense_col based on the nonzero
        // pattern of L(:,ks), and additionally iterate across all nonzeros in
        // Lk_dense_col to find if any entry in Lk_dense_col is untouched, then
        // perform the following
        // L(i,k) = L(i,k)*U(k, Q[ks])/L(P(k),k)
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // pending_scale = sd[k-1]/sd[ks-1]
        if (k > 0)
        {
            SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k-1]));
        }
        else
        {
            SPEX_CHECK(SPEX_mpq_set_ui(pending_scale, 1, 1));
        }
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, sd[ks-1]));
        // remove common factor in mpq_den and mpq_num
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

        // pending_scale *= S(1,ks)
        SPEX_CHECK(SPEX_mpq_mul(pending_scale, pending_scale, SL(ks)));

        // U(k,Q[ks]) will be the new k-th pivot, so we update sd[k] and insert
        // it as the first nonzero in L->v[k]
        // sd[k] = U(k,Q(ks))*S(2,k)
        SPEX_CHECK(SPEX_mpz_divexact(sd[k], Uk_dense_row->x[Qks],
                                     SPEX_MPQ_DEN(SU(k))));
        SPEX_CHECK(SPEX_mpz_mul(sd[k], sd[k], SPEX_MPQ_NUM(SU(k))));
        // update L(P(k),k) = sd[k], which will be the first entry in the list
        SPEX_CHECK(SPEX_mpz_set(L->v[k]->x[0], sd[k]));
        L->v[k]->i[0] = Pk;
        Lk_nz = 1; // Lk_nz is the # of nnz in L->v[k]

        // perform backtracking for each nonzero in col ks of L and store
        // results in Lk_dense_col
        // NOTE: this will cause fillin in the k(th) column of L
        // make sure L->v[k] has enough space to hold all the nnz in L(:,ks) and
        // U(k,ks) that will be inserted
        if (L->v[k]->nzmax < L->v[ks]->nz+1)
        {
            SPEX_CHECK(SPEX_vector_realloc(L->v[k],
                SPEX_MIN(n, SPEX_MAX(2*(L->v[k]->nzmax), L->v[ks]->nz+1)),
                NULL));
        }

        // update entries in L(:,k) for each nnz in L(:,ks)
        // Lk_untouched is the number of entries in Lk_dense_col that have not
        // been used to backtrack yet
        int64_t Lk_untouched = Lk_dense_col->nz;
        for (pks = 0; pks < L->v[ks]->nz; pks++)
        {
            // row index in column ks of L
            cks = L->v[ks]->i[pks];

            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            // L(cks,k) = L(cks,k)*U(k,Q[ks])/L(P[k],k)+L(cks,ks)*pending_scale
            //          = L(cks,k)*  sd[k]   /L(P[k],k)+L(cks,ks)*pending_scale
            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[ks]->x[pks]));
            if (sgn != 0)
            {
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[cks]));
                if (sgn != 0)
                {
#ifndef SPEX_DEBUG
                    // L(cks,k) = floor(L(cks,k)*sd[k]/L(P[k],k))
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[cks],
                        Lk_dense_col->x[cks], sd[k]));
                    SPEX_CHECK(SPEX_mpz_fdiv_q(Lk_dense_col->x[cks],
                        Lk_dense_col->x[cks], Lk_dense_col->x[Pk]));

                    // tmpz = ceil(L(cks,ks)*pending_scale)
                    SPEX_CHECK(SPEX_mpz_mul(tmpz, L->v[ks]->x[pks],
                        SPEX_MPQ_NUM(pending_scale)));
                    SPEX_CHECK(SPEX_mpz_cdiv_q(tmpz, tmpz,
                        SPEX_MPQ_DEN(pending_scale)));
#else
                    // L(cks,k) = floor(L(cks,k)*sd[k]/L(P[k],k))
                    mpq_t r1, r2; mpq_init(r1);mpq_init(r2);
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[cks],
                        Lk_dense_col->x[cks], sd[k]));
                    mpz_fdiv_qr(Lk_dense_col->x[cks],
                        SPEX_MPQ_NUM(r1),
                        Lk_dense_col->x[cks], Lk_dense_col->x[Pk]);
                    mpq_set_den(r1,Lk_dense_col->x[Pk]);
                    mpq_canonicalize(r1);

                    // tmpz = ceil(L(cks,ks)*pending_scale)
                    SPEX_CHECK(SPEX_mpz_mul(tmpz, L->v[ks]->x[pks],
                        SPEX_MPQ_NUM(pending_scale)));
                    mpz_cdiv_qr(tmpz, SPEX_MPQ_NUM(r2), tmpz,
                        SPEX_MPQ_DEN(pending_scale));
                    mpq_set_den(r2, SPEX_MPQ_DEN(pending_scale));
                    mpq_canonicalize(r2);
                    mpq_neg(r2,r2);
                    if (mpq_cmp(r1,r2) != 0)
                    {
        printf("file %s line %d\n",__FILE__,__LINE__);
                        SPEX_CHECK(SPEX_gmp_printf("%Qd\n%Qd\n",r1,r2));
                        mpq_clear(r1);
                        mpq_clear(r2);
                        SPEX_FREE_ALL;
                        return SPEX_PANIC;
                    }
                    mpq_clear(r1);
                    mpq_clear(r2);
#endif
                    // L(cks,k) = L(cks,k)+tmpz
                    SPEX_CHECK(SPEX_mpz_add(Lk_dense_col->x[cks],
                        Lk_dense_col->x[cks], tmpz));
                    Lk_untouched--;
                }
                else  // fillin 
                {
                    // L(cks,k) = L(cks,ks)*pending_scale
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[cks],
                        L->v[ks]->x[pks], SPEX_MPQ_DEN(pending_scale)));
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[cks],
                        Lk_dense_col->x[cks], SPEX_MPQ_NUM(pending_scale)));
                }
            }
            else
            {
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[cks]));
                if (sgn != 0)
                {
                    // L(cks,k) = L(cks,k)*sd[k]/L(P[k],k)
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[cks],
                        Lk_dense_col->x[cks], sd[k]));
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[cks],
                        Lk_dense_col->x[cks], Lk_dense_col->x[Pk]));
                    Lk_untouched--;
                }
            }

            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            // move nonzero entry from Lk_dense_col to L->v[k]->x
            // NOTE: explicit zero due to exact cancellation in backtracking
            //       is removed
            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[cks]));
            if (sgn != 0)
            {
                L->v[k]->i[Lk_nz] = cks;
                SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[Lk_nz],
                                         Lk_dense_col->x[cks]));
                // reset L(cks,k)=0 in the scattered vector
                SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col->x[cks], 0));
                Lk_nz++;
            }
        }

        // continue the backtracking process in case explicit zeros are removed
        // in column ks, which means certain entries in column k remain
        // untouched.
        // check if there is any other nnz left untouched other than L(P[k],k)
        if (Lk_untouched > 1)
        {
            // Lk_untouched is the max # of nnz that would be inserted
            if (L->v[k]->nzmax < Lk_nz+Lk_untouched)
            {
                SPEX_CHECK(SPEX_vector_realloc(L->v[k],
                       SPEX_MIN(n, SPEX_MAX(2*(L->v[k]->nzmax),
                                            Lk_nz+Lk_untouched)), NULL));
            }
            for (pk = 0; pk < Lk_dense_col->nz; pk++)
            {
                // row index in scattered form of column k of L
                ck = Lk_dense_col->i[pk];

                // skip L(P[k],k)
                if (ck == Pk)      {continue;}

                // remove explicit 0s
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[ck]));
                if (sgn == 0)
                {
                    continue;
                }

                // L(ck,k) = L(ck,k)*sd[k]/L(P[k],k)
                SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ck],
                                        Lk_dense_col->x[ck], sd[k]));
                SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ck],
                                             Lk_dense_col->x[ck],
                                             Lk_dense_col->x[Pk]));
                L->v[k]->i[Lk_nz] = ck;
                SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[Lk_nz],
                                         Lk_dense_col->x[ck]));
                // reset L(ck,k)=0 in the scattered vector
                SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col->x[ck], 0));
                Lk_nz++;
            }
        }
        // update number of nnz
        L->v[k]->nz = Lk_nz;
        // For L(P[k],k), we don't need to place it as U(k,Q(ks)) since
        // k-th col will be deleted. Therefore, just reset it to 0.
        SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col->x[Pk], 0));
    }
    else  // U(ks+1:k-1,Q(ks)) contains nnz(s)
    {
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // construct column k of L based on Lk_dense_col
        // explicit zeros are not removed/skipped, since they should be rarely
        // found here (appear only when exact cancellation happens) and also
        // this vector will be updated right after the backtracking process.
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        if (Lk_dense_col->nz > L->v[k]->nzmax)
        {
            SPEX_CHECK(SPEX_vector_realloc(L->v[k],
                SPEX_MIN(n, SPEX_MAX(2*(L->v[k]->nzmax), Lk_dense_col->nz)),
                NULL));
        }

        // put pivot as the first entries in L->v[k]->x
        SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[0], Lk_dense_col->x[Pk]));
        SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col->x[Pk], 0));
        L->v[k]->i[0] = Pk;
        Lk_nz = 1;
        for (pk = 0; pk < Lk_dense_col->nz; pk++)
        {
            ck = Lk_dense_col->i[pk];
            // pivot has been inserted
            if (ck == Pk)            {continue;}
            // swap the entries in the Lk_dense_col and L->v[k]->x
            SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[Lk_nz], Lk_dense_col->x[ck]));
            SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col->x[ck], 0));
            L->v[k]->i[Lk_nz] = ck;
            Lk_nz++;
        }
        L->v[k]->nz = Lk_nz;
        // Lk_dense_col->nz = 0;

        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // initialize history vector h and copy L->v[ks]->x to Lk_dense_col
        // with scale applied. Explicit zero(s) are kept if exist
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // check if S(1,ks) == 1
        SPEX_CHECK(SPEX_mpq_equal(&sgn, SL(ks), one));
        // first entry is the pivot, so we just update it with sd[ks]
        ASSERT(L->v[ks]->i[0] == Pks);
        SPEX_CHECK(SPEX_mpz_set(Lk_dense_col->x[Pks], sd[ks]));
        Lk_dense_col->i[0] = Pks;
        h[Pks] = SPEX_FLIP(ks-1);
        // continue copy the remaining entries
        for (pks = 1; pks < L->v[ks]->nz; pks++) 
        {
            cks = L->v[ks]->i[pks];
            // apply scale to each entry. This must be done to keep these
            // entries in the same scale (with no pending scale) as other
            // entries from U(k:ks-1,Q[ks]), which will have any pending
            // applied as well
            if (sgn == 0) // S(1,ks) != 1
            {
                SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[cks],
                                             L->v[ks]->x[pks],
                                             SPEX_MPQ_DEN(SL(ks))));
                SPEX_CHECK(SPEX_mpz_mul     (Lk_dense_col->x[cks],
                                             Lk_dense_col->x[cks],
                                             SPEX_MPQ_NUM(SL(ks))));
            }
            else
            {
                SPEX_CHECK(SPEX_mpz_set(Lk_dense_col->x[cks],L->v[ks]->x[pks]));
            }
            // set the column index
            Lk_dense_col->i[pks] = cks;
            // Assuming that zero resulted from exact cancellation is kept in
            // the nnz pattern, then performing backtracking will not introduce
            // new fillin. That is, for any i such that U(i,Q[ks])!=0, the nnz
            // pattern of L(i+1:n, i) is a subset of the nnz pattern of
            // Lk_dense_col(i+1:n). Therefore, the final nnz pattern can be
            // found by L->v[ks]->i and column-wise nnz pattern of U(:,ks).
            // We can just initialize the history vector as:
            //
            // h[cks] = ks-1;
            //
            // However, when explicit zero(s) are always eleminated, the
            // following initialization should be used instead:
            //
            // h[cks] = SPEX_FLIP(ks-1); 
            //
            // With such initialization, entry with h > -1 is clearly not in
            // nnz pattern and any entry in the nnz pattern with h = -1 must be
            // nonzero. In all, any explicit zero with h >= -1 must not be in
            // the nnz pattern.  In this way, we can determine if a zero entry
            // in Lk_dense_col is in the nnz pattern.
            h[cks] = SPEX_FLIP(ks-1);
        }
        Lk_nz = L->v[ks]->nz; // # of nnz in Lk_dense_col

        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // backtrack column ks of L and U for each nnz in U(k:ks-1,Q(ks))
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        pks = Uc_ks_nz-1;     // 2nd last entry in the Q[ks]-th col of U
        ASSERT(pks >= 0);
        i = Uci_ks[pks];      // row index
        ASSERT(i > k);
        while (i >= k)
        {
            // Uiks = U(i,Q(ks))*S(2,k)*S(3,k)
            if (i == k)
            {
                // store this pointer for later use
                Ucp_k_ks = pks;

                // Uiks = U(k,Q(ks))*S(2,k)
                SPEX_CHECK(SPEX_mpz_divexact(Uiks, Uk_dense_row->x[Qks],
                                             SPEX_MPQ_DEN(SU(k))));
                SPEX_CHECK(SPEX_mpz_mul(Uiks, Uiks, SPEX_MPQ_NUM(SU(k))));
            }
            else
            {
                ASSERT(U->v[i]->i[Ucx_ks[pks]] == Qks);
                // skip if U(i, Q[ks]) turns out to be explicit zero
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[i]->x[Ucx_ks[pks]]));
                if (sgn == 0)
                {
                    // get next nnz
                    pks --;
                    // the last nonzero must be U(k,Q(ks)), if this is not
                    // in the nonzero pattern found, then it is a new nnz
                    // caused by fillin. Otherwise, cppu won't be used
                    // according to the heuristics in SPEX_LUU
                    if (pks < 0)  {i = k;}
                    else          {i = SPEX_MAX(k, Uci_ks[pks]);}
                    continue;
                }

                // Uiks = U(i,Q(ks))*S(2,i)
                SPEX_CHECK(SPEX_mpz_divexact(Uiks, U->v[i]->x[Ucx_ks[pks]],
                                             SPEX_MPQ_DEN(SU(i))));
                SPEX_CHECK(SPEX_mpz_mul(Uiks, Uiks, SPEX_MPQ_NUM(SU(i))));
            }

            // r = (S(1,i) == 1)
            SPEX_CHECK(SPEX_mpq_equal(&r, SL(i), one));

            // the pivot L(P[i],i) is the first entry in the L->v[i], which
            // should be skipped, so we start with pi = 1
            Pi = P[i];
            for (pi = 1; pi < L->v[i]->nz; pi++)
            {
                // row index of entry in column i of L
                ci = L->v[i]->i[pi];
                ASSERT(ci != Pi);// L(P[i],i) has been excluded

                // Lk_dense_col = (Lk_dense_col*sd(i-1) + L(:,i)*Uiks)/sd(i)
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[ci]));
                if (sgn != 0)
                {
                    ASSERT(h[ci] < -1);
                    if (i > 0)
                    {
                        // Lk_dense_col = Lk_dense_col*sd(i-1)
                        SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ci],
                                                Lk_dense_col->x[ci], sd[i-1]));
                    }
                    int64_t real_hci = SPEX_FLIP(h[ci]);
                    if (i < real_hci || r == 0)
                    {
                        // Lk_dense_col = Lk_dense_col/sd(h[ci])
                        //                + L(ci,i)*Uiks/L(P(i),i)
                        // use L(P(i),i) instead of sd[i] to avoid scaling
                        // for L(:,i)
                        SPEX_CHECK(SPEX_mpz_mul(tmpz, L->v[i]->x[pi], Uiks));
                        SPEX_CHECK(SPEX_mpz_cdiv_q(tmpz, tmpz,
                                   L->v[i]->x[0]));
                        SPEX_CHECK(SPEX_mpz_fdiv_q(Lk_dense_col->x[ci],
                                   Lk_dense_col->x[ci], sd[real_hci]));

                        SPEX_CHECK(SPEX_mpz_add(Lk_dense_col->x[ci],
                                   Lk_dense_col->x[ci], tmpz));
                    }
                    else
                    {
                        // Lk_dense_col = (Lk_dense_col + L(ci,i)*Uiks)/sd(i)
                        SPEX_CHECK(SPEX_mpz_addmul(Lk_dense_col->x[ci],
                                                   L->v[i]->x[pi], Uiks));
                        SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ci],
                                                   Lk_dense_col->x[ci], sd[i]));
                    }
                }
                else
                {
                    if (h[ci] >= -1)
                    {
                        // this is a fill-in
                        Lk_dense_col->i[Lk_nz] = ci;
                        Lk_nz++;
                    }
                    // Lk_dense_col =  L(ci,i)*Uiks/L(P(i),i)
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ci],
                                            L->v[i]->x[pi], Uiks));
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ci],
                                                 Lk_dense_col->x[ci],
                                                 L->v[i]->x[0]));
                }

                // update h[ci]
                h[ci] = SPEX_FLIP(i-1);
            }

            if (i > k)
            {
                // move Uiks (scaled U(i,Q(ks)) to Lk_dense_col
                SPEX_CHECK(SPEX_mpz_swap(Lk_dense_col->x[Pi], Uiks));
                Lk_dense_col->i[Lk_nz] = Pi;
                Lk_nz++;
                // update corresponding entry in the history vector
                h[Pi] = SPEX_FLIP(i-1);
            }
            else //i == k, which is the last loop
            {
                // Uiks will be the new pivot of column k, so we just place it
                // as the first entry in L->v[k]. However, we postpone updating
                // sd[k], since the original value of sd[k] could be used
                SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[0], Uiks));
                L->v[k]->i[0] = Pi;
                break;
            }

            // get next nnz
            pks --;
            // the last nonzero must be U(k,Q(ks)), if this is not
            // in the nonzero pattern found, then it is a new nnz
            // caused by fillin. Otherwise, cppu won't be used
            // according to the heuristics
            if (pks < 0)  {i = k;}
            else          {i = SPEX_MAX(k, Uci_ks[pks]);}
        }
        Lk_dense_col->nz = Lk_nz;
        GOTCHA;

        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // 1. Iterate across all nnz in Lk_dense_col, perform history update if
        //    needed, then move all nonzero entry from Lk_dense_col to
        //    L->v[k]->x
        // NOTE: explicit zero due to exact cancellation in backtracking
        //       is removed.
        // 2. Swap values from L->v[ks]->x and Lk_dense_col
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // reallocate the nonzero pattern if needed
        if (L->v[k]->nzmax < Lk_dense_col->nz+1)
        {
            SPEX_CHECK(SPEX_vector_realloc(L->v[k],
                   SPEX_MIN(n, SPEX_MAX(2*(L->v[k]->nzmax),
                                        Lk_dense_col->nz+1)), NULL));
        }
        Lk_nz = 1; // the pivot has been inserted to L->v[k]
        for (pks = 0; pks < Lk_dense_col->nz; pks++) 
        {
            cks = Lk_dense_col->i[pks];
            h[cks] = SPEX_FLIP(h[cks]);
            ASSERT(h[cks] >= -1);
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[cks]));
            if (sgn != 0)
            {
                // check if need to perform history update
                if (h[cks] != k-1)
                {
                    if (k > 0)
                    {
                        SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[cks],
                                   Lk_dense_col->x[cks], sd[k-1]));
                    }
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[cks],
                               Lk_dense_col->x[cks], sd[h[cks]]));
                }
                L->v[k]->i[Lk_nz] = cks;
                SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[Lk_nz],
                                         Lk_dense_col->x[cks]));
                SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col->x[cks], 0));
                Lk_nz++;
            }
        }
        // update number of nnz in column k of L
        L->v[k]->nz = Lk_nz;
        // update sd[k] = L->v[k]->x[0]
        SPEX_CHECK(SPEX_mpz_set(sd[k], L->v[k]->x[0]));
    }
    // update S(1,k) and S(2,k) as 1, since all entry in L(:,k) are scaled
    SPEX_CHECK(SPEX_mpq_set(SL(k), one));

    //-------------------------------------------------------------------------
    // swap values from L->v[ks]->x and Lk_dense_col
    //-------------------------------------------------------------------------
    pk = 0;
    for (pks = 0; pks < L->v[ks]->nz; pks++)
    {
        cks = L->v[ks]->i[pks];
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[ks]->x[pks]));
        if (sgn != 0 && P_inv[cks] < *inext && P_inv[cks] > ks)
        {
            *inext = P_inv[cks];
        }
        SPEX_CHECK(SPEX_mpz_swap(Lk_dense_col->x[cks], L->v[ks]->x[pks]));
        Lk_dense_col->i[pk] = cks;
        pk++;
    }
    Lk_dense_col->nz = pk;

    //-------------------------------------------------------------------------
    // remove explicit zeros from Uk_dense_row
    //-------------------------------------------------------------------------
    int64_t Uk_nz;
#ifdef SPEX_DEBUG
    Uk_nz = Uk_dense_row->nz;
    for (pk = 0; pk < Uk_nz;)
    {
        j = Uk_dense_row->i[pk];
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row->x[j]));
        ASSERT(sgn != 0);
        // There is no need to remove explicit zeros here, since every time when
        // entries are swapped from U->v[k] to Uk_dense_row, explicit zeros are
        // removed. In addition, cppu is never called after dppu2, while
        // Uk_dense_row is not modified in dppu1, and updated at the end of
        // every call to cppu.
        /*if (sgn == 0)
        {
            // remove explicit zeros
            Uk_nz--;
            Uk_dense_row->i[pk] = Uk_dense_row->i[Uk_nz];
            continue;
        }*/

        // there should not be any nnz in U(k,Q(k+1:ks-1))
        ASSERT (Q_inv[j] == k || Q_inv[j] >= ks);
        pk++;
    }
    Uk_dense_row->nz = Uk_nz;
#endif

    //-------------------------------------------------------------------------
    // RwSOP. This could cause fill-ins!
    //-------------------------------------------------------------------------
    int64_t p_Uiks;  // pointer for U(i,Q(ks))
    // pending_scale = U(k, Q(ks))/ U(k, Q(k))
    SPEX_CHECK(SPEX_mpq_set_z(pending_scale, Uk_dense_row->x[Qks]));
    SPEX_CHECK(SPEX_mpq_set_den(pending_scale, Uk_dense_row->x[Qk]));
    SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

    // iterate for all nnz in U(k+1:ks,Q(ks))
    for (pks = Ucp_k_ks+1; pks <= Uc_ks_nz; pks++)
    {
        // U(ks, Q[ks]) is not in the nnz pattern provided by Uci_ks, Ucx_ks
        if (pks == Uc_ks_nz)
        {
            i = ks;
            p_Uiks = 0;
        }
        else
        {
            i = Uci_ks[pks];       // row index
            p_Uiks = Ucx_ks[pks];  // pointer for U(i,Q(ks))
        }
        int64_t Ui_nz = U->v[i]->nz;

        // skip scaling for frames between iterations
        // REMARK: U(k,Q(ks)) may not be in the nnz pattern initially
        int64_t i1 = (pks == 0) ? k+1 : SPEX_MAX(k, Uci_ks[pks-1])+1;
        for (; i1 < i; i1++)
        {
            SPEX_CHECK(SPEX_mpq_mul(SL(i1), SL(i1), pending_scale));
            SPEX_CHECK(SPEX_mpq_mul(SU(i1), SU(i1), pending_scale));
            SPEX_CHECK(SPEX_mpz_divexact(sd[i1], sd[i1],
                                    SPEX_MPQ_DEN(pending_scale)));
            SPEX_CHECK(SPEX_mpz_mul(sd[i1], sd[i1],
                                    SPEX_MPQ_NUM(pending_scale)));
        }

        // simply scale up whole frame i if U(i, Q[ks]) == 0
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[i]->x[p_Uiks]));
        if (sgn == 0)
        {
            // remove U(i, Q[ks])
            Ui_nz--;
            U->v[i]->i[p_Uiks] = U->v[i]->i[Ui_nz];
            SPEX_CHECK(SPEX_mpz_swap(U->v[i]->x[p_Uiks], U->v[i]->x[Ui_nz]));
            U->v[i]->nz = Ui_nz;

            // skip scaling if U(i, Q(ks)) == 0
            SPEX_CHECK(SPEX_mpq_mul(SL(i), SL(i), pending_scale));
            SPEX_CHECK(SPEX_mpq_mul(SU(i), SU(i), pending_scale));
            SPEX_CHECK(SPEX_mpz_divexact(sd[i], sd[i],
                                    SPEX_MPQ_DEN(pending_scale)));
            SPEX_CHECK(SPEX_mpz_mul(sd[i], sd[i],
                                    SPEX_MPQ_NUM(pending_scale)));
            continue;
        }

        // This is used to determine the number of additional fillins that will
        // be introduced to row i of U after RwSOP. Initializing with
        // Uk_dense_row->nz-2 is because that U(k,Q[k]) and U(k,Q[ks]) won't be
        // used. This number decreases for every other nonzero entries in
        // U(k,:) used to perform RwSOP update for the corresponding entry in
        // U(i,:). 
        int64_t num_of_fillin = Uk_dense_row->nz-2;

        // the updates for all entries in row i of U (except U(i,Q[ks])) share
        // one common feature, that they all need to divide U(k,Q[k]).
        // Therefore, we set
        // S(2, i) = S(2, i)/ U(k,Q(k))
        SPEX_CHECK(SPEX_mpz_mul(SPEX_MPQ_DEN(SU(i)),
                                SPEX_MPQ_DEN(SU(i)),
                                Uk_dense_row->x[Qk]));
        SPEX_CHECK(SPEX_mpq_canonicalize(SU(i)));

        // store the value for U(i,Q[ks])
        ASSERT(U->v[i]->i[p_Uiks] == Qks);
        SPEX_CHECK(SPEX_mpz_set_ui(Uiks, 0));
        SPEX_CHECK(SPEX_mpz_swap(Uiks, U->v[i]->x[p_Uiks]));

        // update row i of U
        // for U(ci,i) with Q_inv(ci) < ks, multiply U(k,Q(ks))
        // for U(ci,i) with Q_inv(ci) > ks but U(k,ci) == 0, multiply U(k,Q(ks))
        // for U(ci,i) with Q_inv(ci) > ks but U(k,ci) != 0, perform 
        //     U(i,ci) = U(i,ci)*U(k,Q(ks)) - U(i,Q(ks))*U(k,ci)
        // then divide all U(ci,i) by the denominator of S(2,i)
        // update U(i,Q[ks]) after iteration, which only needs flipping sign
        for (pi = 0; pi < Ui_nz; pi++)
        {
            ci = U->v[i]->i[pi];
            int64_t real_ci = Q_inv[ci];
            if (real_ci == ks)
            {
                // handle U(i, Q[ks]) after iteration
                h[ci] = -2; // or SPEX_MARK(P, ci);
                continue;
            }
            else if (real_ci > ks)
            {
                // if U(k,ci) is zero then U(i,ci) needs only scaling
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row->x[ci]));
            }
            else
            {
                sgn = 0;
            }

            // multiply U(k, Q[ks])
            SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[pi], U->v[i]->x[pi],
                                    Uk_dense_row->x[Qks]));
            if (sgn != 0)
            {
                // perform RwSOP to U(i,Q(ks+1:n+1))
                // U(i,ci)= U(i,ci)*U(k,Q(ks)) - U(i,Q(ks))*U(k,ci)
                SPEX_CHECK(SPEX_mpz_submul(U->v[i]->x[pi], Uiks,
                                        Uk_dense_row->x[ci]));
                num_of_fillin --;
                h[ci] = -2; // or SPEX_MARK(P, ci);
            }
            SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[pi], U->v[i]->x[pi],
                                    SPEX_MPQ_DEN(SU(i))));
        }

        // flip sign in advance: Uiks = -Uiks
        SPEX_CHECK(SPEX_mpz_neg(Uiks, Uiks));
        // finish RwSOP by checking if there is fillin that should be added
        if (num_of_fillin > 0)
        {
            // allocate additional space if needed
            if (U->v[i]->nzmax < Ui_nz+num_of_fillin)
            {
                SPEX_CHECK(SPEX_vector_realloc(U->v[i],
                       SPEX_MIN(n, SPEX_MAX(2*(U->v[i]->nzmax),
                                            Ui_nz+num_of_fillin)), NULL));
            }
        }
        // add FILLIN and restore P
        for (pk = 0; pk < Uk_dense_row->nz; pk++)
        {
            ck = Uk_dense_row->i[pk];
            if (ck == Qk) { continue; }
            else if (h[ck] != -2) // or (!SPEX_MARKED(P, ck))
            {
                U->v[i]->i[Ui_nz] = ck;
                // U(i,ck)= -U(i,Q(ks))*U(k,ck)
                //        = Uiks       *U(k,ck)
                SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[Ui_nz],
                        Uiks      , Uk_dense_row->x[ck]));

                SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[Ui_nz],
                    U->v[i]->x[Ui_nz],SPEX_MPQ_DEN(SU(i))));
                Ui_nz++;
            }
            else
            {
                h[ck] = -1; // or SPEX_MARK(P, ck);
            }
        }

        // CPPU shortcuts
        if (i != ks)
        {
            // Mathematically, we should update U(i, Q(ks)). However, since
            // column ks of U will become column k after permutation, which
            // will be deleted when finished, we will instead delete
            // U(i,Q(ks)).
            Ui_nz--;
            U->v[i]->i[p_Uiks] = U->v[i]->i[Ui_nz];
            SPEX_CHECK(SPEX_mpz_swap(U->v[i]->x[p_Uiks], U->v[i]->x[Ui_nz]));
            // skip scaling for L(:,i) by setting S(1,i) = S(1,i)*pending_scale
            SPEX_CHECK(SPEX_mpq_mul(SL(i), SL(i), pending_scale));
            // sd[i] = sd[i]*pending_scale
            SPEX_CHECK(SPEX_mpz_divexact(sd[i], sd[i],
                                    SPEX_MPQ_DEN(pending_scale)));
            SPEX_CHECK(SPEX_mpz_mul(sd[i], sd[i],
                                    SPEX_MPQ_NUM(pending_scale)));
        }
        else
        {
            // Since U(ks,Q[ks]) should not change other than flipping sign, and
            // we have let S(2,ks) *= 1/U(k,Q[k]), we need to update it so that
            // it still shares common pending factor S(2,ks) as other entries in
            // row ks.
            // U(ks, Q[ks]) = U(ks, Q[ks]) * U(k,Q(k))
            SPEX_CHECK(SPEX_mpz_mul(Uiks, Uiks, Uk_dense_row->x[Qk]));
            SPEX_CHECK(SPEX_mpz_divexact(Uiks, Uiks, SPEX_MPQ_DEN(SU(i))));

            // move Uiks into row i of U. noted that sign has been flipped!
            SPEX_CHECK(SPEX_mpz_swap(U->v[i]->x[p_Uiks], Uiks));
            SPEX_CHECK(SPEX_mpq_neg(SL(i), SL(i)));
            SPEX_CHECK(SPEX_mpz_neg(sd[i], sd[i]));
        }
        U->v[i]->nz = Ui_nz;
        SPEX_CHECK(SPEX_mpz_set_ui(SPEX_MPQ_DEN(SU(i)), 1));
    }
    GOTCHA;

    //-------------------------------------------------------------------------
    // copy nnz from Uk_dense_row to U(k,:)
    // and copy U(ks,:) to Uk_dense_row
    //-------------------------------------------------------------------------
    // Except the case when swapping with the inserted column (ks == n), cppu
    // can never be used after dppu2 is called. And we know that dppu1 does not
    // modify Uk_dense_row, while at the end of every cppu call, Uk_dense_row
    // is updated with U->v[ks] (which will be U->v[k] after permutation).
    // Therefore, Uk_dense_row always holds the same entries as they were in
    // U->v[k], and the following will always hold.
    ASSERT(Uk_dense_row->nz-1 <= U->v[k]->nzmax);
    /*if (Uk_dense_row->nz-1 > U->v[k]->nzmax)
    {
        SPEX_CHECK(SPEX_vector_realloc(U->v[k], Uk_dense_row->nz-1));
    }*/

    // place the new pivot U(k,Q[ks]) as the first entry
    SPEX_CHECK(SPEX_mpz_swap(U->v[k]->x[0], Uk_dense_row->x[Qks]));
    U->v[k]->i[0] = Qks;
    Uk_nz = 1;
    for (pk = 0; pk < Uk_dense_row->nz; pk++)
    {
        j = Uk_dense_row->i[pk];
        if (j == Qk || j == Qks)
        {
            // no need to copy U(k,Q[k]) since column k will be deleted anyway,
            // while U(k, Q[ks]) has been handled
            SPEX_CHECK(SPEX_mpz_set_ui(Uk_dense_row->x[j], 0));
            continue;
        }

        SPEX_CHECK(SPEX_mpz_swap(U->v[k]->x[Uk_nz], Uk_dense_row->x[j]));
#ifdef SPEX_DEBUG
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row->x[j]));
        ASSERT(sgn == 0);
#endif
        U->v[k]->i[Uk_nz] = j;
        Uk_nz++;
    }
    U->v[k]->nz = Uk_nz;
    GOTCHA;

    // no need to do this when ks>=n-1
    if (ks < n-1)
    {
        Uk_nz = 0;
        // construct a scattered vector for U->v[ks]
        for (pks = 0; pks < U->v[ks]->nz; pks++)
        {
            j = U->v[ks]->i[pks];
            int64_t real_j = Q_inv[j];
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[ks]->x[pks]));
            if (sgn == 0) {continue;}

            if (real_j == ks)
            {
                // this entry should be considered as the IPGE update of the
                // entry in the Q[k]-th column
                j = Qk;
            }
            else if (real_j < *jnext)
            {
                // update the index for the next nnz in current row
                *jnext = real_j;
            }
            Uk_dense_row->i[Uk_nz] = j;
            SPEX_CHECK(SPEX_mpz_swap(Uk_dense_row->x[j], U->v[ks]->x[pks]));
            Uk_nz++;
        }
        Uk_dense_row->nz = Uk_nz;
    }
    GOTCHA;

    //-------------------------------------------------------------------------
    // update column permutation
    //-------------------------------------------------------------------------
    Q[k] = Qks;          Q[ks] = Qk;
    Q_inv[Qks] = k;  Q_inv[Qk] = ks;
    GOTCHA;

    //-------------------------------------------------------------------------
    // flip sign for columns and rows ks+1 to n and update corresponding sd
    //-------------------------------------------------------------------------
    for (i = ks+1; i < n; i++)
    {
        SPEX_CHECK(SPEX_mpq_neg(SL(i), SL(i)));
        SPEX_CHECK(SPEX_mpq_neg(SU(i), SU(i)));
        SPEX_CHECK(SPEX_mpz_neg(sd[i], sd[i]));
    }
    GOTCHA;

    SPEX_FREE_ALL;
    return SPEX_OK;
}
