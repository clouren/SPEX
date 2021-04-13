//------------------------------------------------------------------------------
// SPEX_CHOLMOD/spex_cppu.c: perform column permutation pivot update
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

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
    SPEX_MPQ_CLEAR(tmpq);            \
    SPEX_MPQ_CLEAR(pending_scale);   \
    SPEX_MPZ_CLEAR(tmpz);

#include "spex_lu_update_internal.h"

SPEX_info spex_cppu
(
    SPEX_mat *L,  // matrix L
    SPEX_mat *U,  // matrix U
    mpq_t *S,        // array of size 3*n that stores pending scales
    mpz_t *d,        // array of size n that stores the unscaled pivot
    mpz_t *sd,       // array of size n that stores the scaled pivot
    spex_scattered_vector *Lk_dense_col,// scattered column k of L
    spex_scattered_vector *Uk_dense_row,// scattered column k of U
    int64_t *inext,  // the index of first off-diag entry in col k of L
    int64_t *jnext,  // the index of first off-diag entry in row k of U
    int64_t *h,      // allocated vector that can be used for history vector.
                     // All entries are maintained to be >= -1
    int64_t *Q,      // column permutation
    int64_t *Q_inv,  // inverse of column permutation
    const int64_t *P,// row permutation
    const int64_t *P_inv,// inverse of row permutation
    int64_t *Ldiag,  // L(P(k),k) can be found as L->v[k]->x[Ldiag[k]]
    const int64_t *Uci,// the row index for col-wise nnz pattern of U
    const int64_t *Ucp,// col pointers for col-wise nnz pattern of U
    const int64_t *Ucx,// the value of k-th entry is found as
                     // U->v[Uci[k]]->x[Ucx[k]]
    const int64_t k, // current column index 0 <= k < n
    const int64_t ks // index of the diagonal to be swapped with, [0,n)
)
{
    printf("using cppu swapping k(%ld) and ks(%ld)\n",k,ks);
    // initialize workspace
    SPEX_info info;
    int sgn, r;
    int64_t pk, ck, pks, cks, pi, ci, i, j, n = U->n;
    // the pointer for U(k,Q(ks)) = U->v[k]->x[Ucx[Ucp_k_ks]]
    int64_t Ucp_k_ks;
    *inext = n;
    *jnext = n;

    mpq_t pending_scale, tmpq, one;
    SPEX_MPQ_SET_NULL(pending_scale); SPEX_MPQ_SET_NULL(tmpq);
    mpz_t Uiks, tmpz; SPEX_MPZ_SET_NULL(Uiks); SPEX_MPZ_SET_NULL(tmpz);
    SPEX_CHECK(SPEX_mpq_init(pending_scale));
    SPEX_CHECK(SPEX_mpq_init(tmpq));
    SPEX_CHECK(SPEX_mpq_init(one));
    SPEX_CHECK(SPEX_mpq_set_ui(one, 1, 1));
    SPEX_CHECK(SPEX_mpz_init(Uiks));
    SPEX_CHECK(SPEX_mpz_init(tmpz));

    if (ks == n)
    {
        printf("ks==n using cppu\n");
        // since the value in Uk_dense_row[Q[k]] will not be used, we use it to
        // hold the original value of sd[k] before swapping column k with
        // column n-1. Then we set sd[k] to d[k]=vk[P[k]]
        SPEX_CHECK(SPEX_mpz_swap(Uk_dense_row->x[Q[k]], sd[k]));
        SPEX_CHECK(SPEX_mpz_set(sd[k], d[k]));

        // get the scale for entries between frames k and n-1
        // pending_scale = sd(k)/Uk_dense_row[Q[k]]
        SPEX_CHECK(SPEX_mpq_set_z(pending_scale, sd[k]));
        SPEX_CHECK(SPEX_mpq_set_den(pending_scale, Uk_dense_row->x[Q[k]]));
        SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

        // if the inserted column is used, we won't need to perform
        // backtracking. Instead, we just need to perform RwSOP.
        //
        // If U(k, Q[n-1]) == 0, RwSOP is simplified as pure scaling for all
        // frame from k+1:n-1. Otherwise (which won't happen due to heuristic),
        // we need to first compute the (n-1)-th IPGE iteration for the column
        // k, and use the result to perform RwSOP on column n-1 of U.
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row->x[Q[n-1]]));
        if (sgn == 0)
        {
            // just need to perform RwSOP by scaling all frame k+1:n-1
            for (j = k+1; j < n; j++)
            {
                // S(3,k+1:n-1) = S(3,k+1:n-1)*pending_scale;
                SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 3, j),
                                        SPEX_LUU_2D(S, 3, j), pending_scale));
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
            // S(2,k) = S(2,k)*S(3,k)
            SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 2, k),
                                    SPEX_LUU_2D(S, 2, k), SPEX_LUU_2D(S, 3, k)));
            // U(k,Q[n-1]) = U(k, Q[n-1])*S(2,k)
            SPEX_CHECK(SPEX_mpz_divexact(Uk_dense_row->x[Q[n-1]],
                      Uk_dense_row->x[Q[n-1]], SPEX_MPQ_DEN(SPEX_LUU_2D(S, 2, k))));
            SPEX_CHECK(SPEX_mpz_mul(Uk_dense_row->x[Q[n-1]],
                      Uk_dense_row->x[Q[n-1]], SPEX_MPQ_NUM(SPEX_LUU_2D(S, 2, k))));
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
            // initialize history vector
            for (pk = 0; pk < Lk_dense_col->nz; pk++)
            {
                ck = Lk_dense_col->i[pk];
                h[ck] = SPEX_FLIP(k-1);
            }
            for (pk = 0; pk < L->v[k]->nz; pk++)
            {
                ck = L->v[k]->i[pk];
                if (ck == P[k])     { continue; }// skip updating L(P(k),k)
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[k]->x[pk]));
                if (sgn == 0)       { continue;   }

                // L(ck,k) = (L(ck, k)*vk(P[k])-L(P[k],k)*vk(ck))/sd[k-1]
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[ck]));
                if (sgn != 0)
                {
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ck],
                                            Lk_dense_col->x[ck], d[k]));
                }
                else if (h[ck] >= -1) //this entry wasn't in the nnz pattern
                {
                    // insert new entry in the nonzero pattern
                    Lk_dense_col->i[Lk_dense_col->nz] = ck;
                    Lk_dense_col->nz++;
                }
                SPEX_CHECK(SPEX_mpz_submul(Lk_dense_col->x[ck],
                                        Lk_dense_col->x[P[k]], L->v[k]->x[pk]));
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
                // S(3,i) = S(3, i)*pending_scale;
                SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 3, i),
                                        SPEX_LUU_2D(S, 3, i), pending_scale));
                // sd[i] = sd[i]*pending_scale
                SPEX_CHECK(SPEX_mpz_divexact(sd[i], sd[i],
                                        SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(sd[i], sd[i],
                                        SPEX_MPQ_NUM(pending_scale)));

                // skip if Lk_dense_col[P[i]] == 0
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[P[i]]));
                if (sgn == 0)       {continue;}

                // perform i-th IPGE update for Lk_dense_col
                SPEX_CHECK(spex_ipge(Lk_dense_col, h,
                    NULL, L->v[i], P, P_inv, (const mpz_t*) sd, d,
                    U->v[i]->x[Ucx[Ucp[Q[i]+1]-1]], SPEX_LUU_2D(S, 1, i),
                    SPEX_LUU_2D(S, 3, i), SPEX_LUU_2D(S, 2, i), i, Ldiag[i]));

                // perform RwSOP for row i with flipped-sign entries in
                // Lk_dense_col. All entries in row i of U must be SCALEUP such
                // that S(:,i)=[S(1,i)*S(3,i);1;1]

                SPEX_CHECK(SPEX_mpq_equal(&r, SPEX_LUU_2D(S, 3, i), one));
                if (r == 0) // S(3,i) != 1
                {
                    // set S(:,i) = [S(1,i)*S(3,i); S(2,i)*S(3,i); 1]
                    SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 1, i),
                                            SPEX_LUU_2D(S, 1, i), SPEX_LUU_2D(S, 3,i)));
                    SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 2, i),
                                            SPEX_LUU_2D(S, 2, i), SPEX_LUU_2D(S, 3,i)));
                    SPEX_CHECK(SPEX_mpq_set_ui(SPEX_LUU_2D(S, 3, i), 1, 1));
                }

                int64_t p = -1; // the pointer to U(i,Q(n-1))
                // iterate all nnz in row i of U
                for (pi = 0; pi < U->v[i]->nz; pi++)
                {
                    ci = U->v[i]->i[pi];
                    if (Q_inv[ci] == n-1)
                    {
                        p = pi;
                    }
                    else if (Q_inv[ci] == i)
                    {
                        // set U(i, Q(i)) = sd[i]
                        SPEX_CHECK(SPEX_mpz_set(U->v[i]->x[pi], sd[i]));
                    }
                    else //if (Q_inv[ci] < n-1)
                    {
                        // apply S(2,i) to U(i,Q(i:n-2))
                        SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[pi],
                            U->v[i]->x[pi], SPEX_MPQ_DEN(SPEX_LUU_2D(S, 2, i))));
                        SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[pi],
                            U->v[i]->x[pi], SPEX_MPQ_NUM(SPEX_LUU_2D(S, 2, i))));
                    }
                }
                // perform RwSOP to U(i,Q(n-1)), POSSIBLE FILLIN
                // sign is changed here due to column swap
                // U(i,ci)= U(i,ci)*S(2,i) +
                //        Lk_dense_col(P[i])*U(k,ci)/Lk_dense_col[P[k]]
                ci = Q[n-1];
                if (p > -1)
                {
                    SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[p],
                               U->v[i]->x[p], SPEX_MPQ_NUM(SPEX_LUU_2D(S, 2, i))));
                    SPEX_CHECK(SPEX_mpz_fdiv_q(U->v[i]->x[p],
                               U->v[i]->x[p], SPEX_MPQ_DEN(SPEX_LUU_2D(S, 2, i))));
                    SPEX_CHECK(SPEX_mpz_mul(tmpz,
                               Lk_dense_col->x[P[i]], Uk_dense_row->x[ci]));
        printf("file %s Line %d\n",__FILE__,__LINE__);
        SPEX_CHECK(SPEX_gmp_printf("%Zd\n",Lk_dense_col->x[P[k]]));
                    SPEX_CHECK(SPEX_mpz_cdiv_q(tmpz, tmpz,
                               Lk_dense_col->x[P[k]]));
        printf("file %s Line %d\n",__FILE__,__LINE__);
                    SPEX_CHECK(SPEX_mpz_add(U->v[i]->x[p],U->v[i]->x[p], tmpz));
                }
                else // U(i,Q(n-1)) was not in the nnz pattern
                {
                    p = U->v[i]->nz;
                    // reallocate the nonzero pattern if needed
                    if (p == U->v[i]->nzmax)
                    {
                        SPEX_CHECK(SPEX_vector_realloc(U->v[i],
                            SPEX_MIN(n, 2*(U->v[i]->nzmax))));
                    }
                    // insert new entry in the nonzero pattern
                    U->v[i]->i[p] = ci;
                    U->v[i]->nz++;

                    SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[p],
                                  Lk_dense_col->x[P[i]], Uk_dense_row->x[ci]));
                    SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[p],
                                  U->v[i]->x[p], Lk_dense_col->x[P[k]]));
                }
                SPEX_CHECK(SPEX_mpq_set_ui(SPEX_LUU_2D(S, 2, i), 1, 1));
                
                if (i == n-1)
                {
                    // make sure U(n-1,Q(n-1)) != 0
                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[n-1]->x[0]));
                    if (sgn == 0)
                    {
                        SPEX_FREE_ALL;
                        return SPEX_SINGULAR;
                    }

                    // L(P(n-1),n-1) = d(n-1) = sd(n-1) = U(n-1, Q(n-1))
                    SPEX_CHECK(SPEX_mpz_set(L->v[n-1]->x[0], U->v[n-1]->x[0]));
                    SPEX_CHECK(SPEX_mpz_set(   d[n-1]      , U->v[n-1]->x[0]));
                    SPEX_CHECK(SPEX_mpz_set(  sd[n-1]      , U->v[n-1]->x[0]));
                    // S(:,n-1) = ones
                    SPEX_CHECK(SPEX_mpq_set(SPEX_LUU_2D(S, 1, n-1), one));
                    SPEX_CHECK(SPEX_mpq_set(SPEX_LUU_2D(S, 2, n-1), one));
                    SPEX_CHECK(SPEX_mpq_set(SPEX_LUU_2D(S, 3, n-1), one));
                }
            }

            // move data from Uk_dense_row, there is only one entry that needs
            // to move, which is U(k,Q[n-1]) and has no pending scale
            SPEX_CHECK(SPEX_mpz_swap(U->v[k]->x[0], Uk_dense_row->x[Q[n-1]]));
            U->v[k]->i[0] = Q[n-1];
            U->v[k]->nz = 1;
        }

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
    // initialized to be 2nd last entry in the Q[ks]-th col of U
    Ucp_k_ks = Ucp[Q[ks]+1]-2;
    // When cppu is called, we know U(k,Q[ks]) must be nnz, while it could be a
    // new fillin added in the update process. If U(k,Q[ks]) exists in the nnz
    // pattern initially, and 2nd last nnz in U(:,Q(ks)) is at row k (i.e.,
    // U(k,Q[ks]) is the 2nd last nnz in column of Q[ks]), then Uci[Ucp_k_ks]
    // == k means that U(k+1:ks-1,Q(ks)) are all zero(s). On the other hand, if
    // U(k,Q[ks]) is a new fillin, then Uci[Ucp_k_ks] < k means that
    // U(k+1:ks-1,Q(ks)) are all zero(s).
    if (Uci[Ucp_k_ks] <= k)
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
        // S(2,k) = S(2,k)*S(3,k)
        SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 2, k), SPEX_LUU_2D(S, 2, k),
                                SPEX_LUU_2D(S, 3, k)));
        // Uiks = U(k,Q(ks))*S(2,k)
        SPEX_CHECK(SPEX_mpz_divexact(Uiks, Uk_dense_row->x[Q[ks]],
                                SPEX_MPQ_DEN(SPEX_LUU_2D(S, 2, k))));
        SPEX_CHECK(SPEX_mpz_mul(Uiks, Uiks,
                                SPEX_MPQ_NUM(SPEX_LUU_2D(S, 2, k))));

        // update d[k] = U(k, Q(ks))
        SPEX_CHECK(SPEX_mpz_set(d[k], Uk_dense_row->x[Q[ks]]));

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

        // pending_scale *= S(1,ks)*S(3,ks)
        SPEX_CHECK(SPEX_mpq_mul(pending_scale, pending_scale, SPEX_LUU_2D(S,1,ks)));
        SPEX_CHECK(SPEX_mpq_mul(pending_scale, pending_scale, SPEX_LUU_2D(S,3,ks)));

        // perform backtracking for each nonzero in col ks of L and store
        // results in Lk_dense_col
        // NOTE: this will cause fillin in the k(th) column of L
        // make sure L->v[k] has enough space to hold all the nnz in L(:,ks) and
        // U(k,ks) that will be inserted
        if (L->v[k]->nzmax < L->v[ks]->nz+1)
        {
            SPEX_CHECK(SPEX_vector_realloc(L->v[k], L->v[ks]->nz+1));
        }

        // update entries in L(:,k) for each nnz in L(:,ks)
        // Lk_nz is the number of nnz currently in L->v[k], and Lk_untouched
        // is the number of entries in Lk_dense_col that have not been used to
        // backtrack yet
        int64_t Lk_nz = 0, Lk_untouched = Lk_dense_col->nz;
        for (pks = 0; pks < L->v[ks]->nz; pks++)
        {
            // row index in column ks of L
            cks = L->v[ks]->i[pks];

            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            // L(cks,k) = L(cks,k)*Uiks/L(P[k],k)+L(cks,ks)*pending_scale
            // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[ks]->x[pks]));
            if (sgn != 0)
            {
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[cks]));
                if (sgn != 0)
                {
                    // L(cks,k) = floor(L(cks,k)*Uiks/L(P[k],k))
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[cks],
                        Lk_dense_col->x[cks], Uiks));
                    SPEX_CHECK(SPEX_mpz_fdiv_q(Lk_dense_col->x[cks],
                        Lk_dense_col->x[cks], Lk_dense_col->x[P[k]]));

                    // tmpz = ceil(L(cks,ks)*pending_scale)
                    SPEX_CHECK(SPEX_mpz_mul(tmpz, L->v[ks]->x[pks],
                        SPEX_MPQ_NUM(pending_scale)));
                    SPEX_CHECK(SPEX_mpz_cdiv_q(tmpz, tmpz,
                        SPEX_MPQ_DEN(pending_scale)));
                    // L(cks,k) = L(cks,k)+tmpz
                    SPEX_CHECK(SPEX_mpz_add(Lk_dense_col->x[cks],
                        Lk_dense_col->x[cks], tmpz));
                    //h[cks] = -2;// mark in history vector
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
                    // L(cks,k) = L(cks,k)*Uiks/L(P[k],k)
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[cks],
                        Lk_dense_col->x[cks], Uiks));
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[cks],
                        Lk_dense_col->x[cks], Lk_dense_col->x[P[k]]));
                    //h[cks] = -2; // mark in history vector
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
                Lk_nz++;
#ifdef SPEX_DEBUG
                // This should become zero after swapping
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[cks]));
                ASSERT(sgn == 0);
#endif
            }
        }

        // continue the backtracking process in case explicit zeros are removed
        // in column ks, which means certain entries in column k remain
        // untouched.
        if (Lk_untouched > 1) // there must be 1 nnz left untouched, L(P[k],k)
        {
            // Lk_untouched is the max # of nnz that would be inserted
            if (L->v[k]->nzmax < Lk_nz+Lk_untouched)
            {
                SPEX_CHECK(SPEX_vector_realloc(L->v[k], Lk_nz+Lk_untouched));
            }
            for (pk = 0; /*Lk_untouched > 1 &&*/ pk < Lk_dense_col->nz; pk++)
            {
                // row index in scattered form of column k of L
                ck = Lk_dense_col->i[pk];

                // skip for L(P[k],k)
                if (ck == P[k]) {continue;}

                // update entries that is not in column ks
                /*if (h[ck] == -2) // updated already
                {
                    h[ck] = -1;
                    continue;
                }*/
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[ck]));
                if (sgn == 0)
                {
                    //Lk_untouched--;
                    continue;
                }

                // L(ck,k) = L(ck,k)*Uiks/L(P[k],k)
                SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ck],
                                        Lk_dense_col->x[ck], Uiks));
                SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ck],
                                             Lk_dense_col->x[ck],
                                             Lk_dense_col->x[P[k]]));
                L->v[k]->i[Lk_nz] = ck;
                SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[Lk_nz],
                                         Lk_dense_col->x[ck]));
                Lk_nz++;
                //Lk_untouched--;
            }
        }
        // update L(P(k),k) and reset Lk_dense_col[P[k]]
        L->v[k]->i[Lk_nz] = P[k];
        SPEX_CHECK(SPEX_mpz_set(sd[k], Uiks));
        SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[Lk_nz], Uiks));
        // no need to put L(P[k],k) to U(k,Q(ks)) since k-th col will be deleted
        SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col->x[P[k]], 0));
        Ldiag[k] = Lk_nz;
        Lk_nz++;
        // update number of nnz
        L->v[k]->nz = Lk_nz;

        // get the pointer for U(k,Q(ks)) = U->v[k]->x[Ucx[Ucp_k_ks]]
        // Ucp_k_ks = Ucp[Q[ks]+1]-2; // 2nd last entry in the Q[ks]-th col of U
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
            SPEX_CHECK(SPEX_vector_realloc(L->v[k], Lk_dense_col->nz));
        }
        for (pk = 0; pk < Lk_dense_col->nz; pk++)
        {
            ck = Lk_dense_col->i[pk];
            if (ck == P[k])
            {
                Ldiag[k] = pk;  // update Ldiag[k]
            }
            // swap the entries in the Lk_dense_col and L->v[k]->x
            SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[pk], Lk_dense_col->x[ck]));
            L->v[k]->i[pk] = ck;
        }
        L->v[k]->nz = Lk_dense_col->nz;
        Lk_dense_col->nz = 0;
#ifdef SPEX_DEBUG
        for(pk = 0; pk < n; pk++)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[pk]));
            ASSERT(sgn == 0);
        }
#endif

        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // initialize history vector h and copy L->v[ks]->x to Lk_dense_col
        // with scale applied. Explicit zero(s) are kept if exist
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // pending_scale = S(1,ks)*S(3,ks)
        SPEX_CHECK(SPEX_mpq_mul(pending_scale,
                                SPEX_LUU_2D(S, 1, ks), SPEX_LUU_2D(S, 3, ks)));
        for (pks = 0; pks < L->v[ks]->nz; pks++) 
        {
            cks = L->v[ks]->i[pks]; 
            // copy value from L->v[ks]->x to Lk_dense_col
            if (cks == P[ks])
            {
                SPEX_CHECK(SPEX_mpz_set(Lk_dense_col->x[cks], sd[ks]));
            }
            else
            {
                // apply scale to each entry
                SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[cks],
                                             L->v[ks]->x[pks],
                                             SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul     (Lk_dense_col->x[cks],
                                             Lk_dense_col->x[cks],
                                             SPEX_MPQ_NUM(pending_scale)));
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
        Lk_dense_col->nz = L->v[ks]->nz;

        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // backtrack column ks of L and U for each nnz in U(k:ks-1,Q(ks))
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        pks = Ucp[Q[ks]+1]-2; // 2nd last entry in the Q[ks]-th col of U
        i = Uci[pks];         // row index
        ASSERT(i > k);
        while (i >= k)
        {
            // Uiks = U(i,Q(ks))*S(2,k)*S(3,k)
            if (i == k)
            {
                // store this pointer for later use
                Ucp_k_ks = pks;

                // update S(2,k)
                SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 2, k), SPEX_LUU_2D(S, 2, k),
                                        SPEX_LUU_2D(S, 3, k)));
                // update d[k]
                SPEX_CHECK(SPEX_mpz_set(d[k], Uk_dense_row->x[Q[ks]]));
                SPEX_CHECK(SPEX_mpz_divexact(Uiks, Uk_dense_row->x[Q[ks]],
                                        SPEX_MPQ_DEN(SPEX_LUU_2D(S, 2, k))));
                SPEX_CHECK(SPEX_mpz_mul(Uiks, Uiks,
                                        SPEX_MPQ_NUM(SPEX_LUU_2D(S, 2, k))));
            }
            else
            {
                // skip if U(i, Q[ks]) turns out to be explicit zero
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[i]->x[Ucx[pks]]));
                if (sgn == 0)
                {
                    // get next nnz
                    pks --;
                    // the last nonzero must be U(k,Q(ks)), if this is not
                    // in the nonzero pattern found, then it is a new nnz
                    // caused by fillin. Otherwise, cppu won't be used
                    // according to the heuristics in SPEX_LUU
                    if (pks < 0)  {i = k;}
                    else          {i = SPEX_MAX(k, Uci[pks]);}
                    continue;
                }

                // Uiks = U(i,Q(ks))*S(2,i)*S(3,i)
                SPEX_CHECK(SPEX_mpq_mul(pending_scale, SPEX_LUU_2D(S, 2, i),
                                        SPEX_LUU_2D(S, 3, i)));
                SPEX_CHECK(SPEX_mpz_divexact(Uiks, U->v[i]->x[Ucx[pks]],
                                        SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(Uiks, Uiks,
                                        SPEX_MPQ_NUM(pending_scale)));
            }

            // r = (S(1,i)*S(3,i) == 1)
            SPEX_CHECK(SPEX_mpq_mul(pending_scale, SPEX_LUU_2D(S, 1, i),
                                    SPEX_LUU_2D(S, 3, i)));
            SPEX_CHECK(SPEX_mpq_equal(&r, pending_scale, one));

            for (pi = 0; pi < L->v[i]->nz; pi++)
            {
                // exclude L(P[i], i)
                if (pi == Ldiag[i])     { continue;}

                // row index of entry in column i of L
                ci = L->v[i]->i[pi];
                ASSERT(ci != P[i]);// L(P[i],i) has been excluded

                // Lk_dense_col = (Lk_dense_col*sd(i-1) + L(:,i)*Uiks)/sd(i)
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[ci]));
                if (sgn != 0)
                {
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
                        //                + L(:,i)*Uiks/L(P(i),i)
                        // use L(P(i),i) instead of sd[i] to avoid scaling
                        // for L(:,i)
                        SPEX_CHECK(SPEX_mpz_mul(tmpz, L->v[i]->x[pi], Uiks));
                        SPEX_CHECK(SPEX_mpz_cdiv_q(tmpz, tmpz,
                                   L->v[i]->x[Ldiag[i]]));

                        SPEX_CHECK(SPEX_mpz_fdiv_q(Lk_dense_col->x[ci],
                                   Lk_dense_col->x[ci], sd[real_hci]));

                        SPEX_CHECK(SPEX_mpz_add(Lk_dense_col->x[ci],
                                   Lk_dense_col->x[ci], tmpz));
                    }
                    else
                    {
                        // Lk_dense_col = (Lk_dense_col + L(:,i)*Uiks)/sd(i)
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
                        Lk_dense_col->i[Lk_dense_col->nz] = ci;
                        Lk_dense_col->nz++;
                    }
                    // Lk_dense_col =  L(:,i)*Uiks/L(P(i),i)
                    SPEX_CHECK(SPEX_mpz_mul(Lk_dense_col->x[ci],
                                            L->v[i]->x[pi], Uiks));
                    SPEX_CHECK(SPEX_mpz_divexact(Lk_dense_col->x[ci],
                                                 Lk_dense_col->x[ci],
                                                 L->v[i]->x[Ldiag[i]]));
                }

                // update h[ci]
                h[ci] = SPEX_FLIP(i-1);
            }

            // move Uiks (scaled U(i,Q(ks)) to Lk_dense_col
            SPEX_CHECK(SPEX_mpz_swap(Lk_dense_col->x[P[i]], Uiks));
            Lk_dense_col->i[Lk_dense_col->nz] = P[i];
            Lk_dense_col->nz++;
            // update corresponding entry in the history vector
            h[P[i]] = SPEX_FLIP(i-1);

            if (i == k) {break;} // k is the last loop

            // get next nnz
            pks --;
            // the last nonzero must be U(k,Q(ks)), if this is not
            // in the nonzero pattern found, then it is a new nnz
            // caused by fillin. Otherwise, cppu won't be used
            // according to the heuristics
            if (pks < 0)  {i = k;}
            else          {i = SPEX_MAX(k, Uci[pks]);}
        }

        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // 1. Iterate across all nnz in Lk_dense_col, perform history update if
        //    needed, then move all nonzero entry from Lk_dense_col to
        //    L->v[k]->x
        // NOTE: explicit zero due to exact cancellation in backtracking
        //       is removed.
        // 2. Swap values from L->v[ks]->x and Lk_dense_col
        // . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
        // reallocate the nonzero pattern if needed
        if (L->v[k]->nzmax < Lk_dense_col->nz)
        {
            SPEX_CHECK(SPEX_vector_realloc(L->v[k], Lk_dense_col->nz));
        }
        pk = 0;
        for (pks = 0; pks < Lk_dense_col->nz; pks++) 
        {
            cks = Lk_dense_col->i[pks]; 
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, Lk_dense_col->x[cks]));
            h[cks] = SPEX_FLIP(h[cks]);
            ASSERT(h[cks] >= -1);
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
                if (cks == P[k]) // find the diagnal of L(P[k],k)
                {
                    SPEX_CHECK(SPEX_mpz_set(sd[k], Lk_dense_col->x[cks]));
                    Ldiag[k] = pk;
                }
                L->v[k]->i[pk] = cks;
                SPEX_CHECK(SPEX_mpz_swap(L->v[k]->x[pk], Lk_dense_col->x[cks]));
                SPEX_CHECK(SPEX_mpz_set_ui(Lk_dense_col->x[cks], 0));
                pk++;
            }
        }
        // update number of nnz in column k of L
        L->v[k]->nz = pk;
    }
    // update S(1,k) and S(2,k) as 1, since all entry in L(:,k) are scaled
    SPEX_CHECK(SPEX_mpq_set(SPEX_LUU_2D(S, 1, k), one));
    SPEX_CHECK(SPEX_mpq_set(SPEX_LUU_2D(S, 3, k), one));

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
    // pending_scale = U(k, Q(ks))/ U(k, Q(k))
    SPEX_CHECK(SPEX_mpq_set_z(pending_scale, Uk_dense_row->x[Q[ks]]));
    SPEX_CHECK(SPEX_mpq_set_den(pending_scale, Uk_dense_row->x[Q[k]]));
    SPEX_CHECK(SPEX_mpq_canonicalize(pending_scale));

    // iterate for all nnz in U(k+1:ks,Q(ks))
    for (pks = Ucp_k_ks+1; pks < Ucp[Q[ks]+1]; pks++)
    {
        // This is used to determine the number of additional fillins that will
        // be introduced to row i of U after RwSOP. Initializing with
        // Uk_dense_row->nz-2 is because that U(k,Q[k]) and U(k,Q[ks]) won't be
        // used. This number decreases for every other nonzero entries in
        // U(k,:) used to perform RwSOP update for the corresponding entry in
        // U(i,:). 
        int64_t num_of_fillin = Uk_dense_row->nz-2;
        // row index
        i = Uci[pks];

        // skip scaling for frames between iterations
        // REMARK: U(k,Q(ks)) may not be in the nnz pattern initially
        int64_t i1 = (pks == 0) ? k+1 : SPEX_MAX(k, Uci[pks-1])+1;
        for (; i1 < i; i1++)
        {
            SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 3, i1),
                                    SPEX_LUU_2D(S, 3, i1), pending_scale));
            SPEX_CHECK(SPEX_mpz_divexact(sd[i1], sd[i1],
                                    SPEX_MPQ_DEN(pending_scale)));
            SPEX_CHECK(SPEX_mpz_mul(sd[i1], sd[i1],
                                    SPEX_MPQ_NUM(pending_scale)));
        }

        // skip scaling if U(i, Q(ks)) == 0
        SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[i]->x[Ucx[pks]]));
        if (sgn == 0)
        {
            SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 3, i),
                                    SPEX_LUU_2D(S, 3, i), pending_scale));
            SPEX_CHECK(SPEX_mpz_divexact(sd[i], sd[i],
                                    SPEX_MPQ_DEN(pending_scale)));
            SPEX_CHECK(SPEX_mpz_mul(sd[i], sd[i],
                                    SPEX_MPQ_NUM(pending_scale)));
            continue;
        }

        for (pk = 0; pk < Uk_dense_row->nz; pk++)
        {
            j = Uk_dense_row->i[pk];
            if (Q_inv[j] != k)
            {
                h[j] = -2;
            }
        }

        // update row i of U
        // for U(ci,i) with Q_inv(ci) < ks, apply pending_scale
        // for U(ci,i) with Q_inv(ci) > ks but U(k,ci) == 0, apply pending_scale
        // for U(ci,i) with Q_inv(ci) < ks but U(k,ci) != 0, perform 
        //     U(i,ci) = (U(i,ci)*U(k,Q(ks)) - U(i,Q(ks))*U(k,ci))/U(k,Q(k))
        // update U(i,Q[ks]) after iteration, which only needs flipping sign
        for (pi = 0; pi < U->v[i]->nz; pi++)
        {
            ci = U->v[i]->i[pi];
            if (Q_inv[ci] == ks)
            {
                // handle U(i, Q[ks]) after iteration
                h[ci] = -1;
                continue;
            }
            else if (Q_inv[ci] > ks)
            {
                // if U(k,ci) is zero then U(i,ci) needs only scaling
                SPEX_CHECK(SPEX_mpz_sgn(&sgn, Uk_dense_row->x[ci]));
            }
            else
            {
                sgn = 0;
            }

            if (sgn == 0)
            {
                // apply pending_scale to U(i,Q(i:ks-1))
                SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[pi], U->v[i]->x[pi],
                                        SPEX_MPQ_DEN(pending_scale)));
                SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[pi], U->v[i]->x[pi],
                                        SPEX_MPQ_NUM(pending_scale)));
            }
            else
            {
                // perform RwSOP to U(i,Q(ks+1:n+1))
                // U(i,ci)= (U(i,ci)*U(k,Q(ks)) - U(i,Q(ks))*U(k,ci))/U(k,Q(k))
                SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[pi], U->v[i]->x[pi],
                                        Uk_dense_row->x[Q[ks]]));
                SPEX_CHECK(SPEX_mpz_submul(U->v[i]->x[pi], U->v[i]->x[Ucx[pks]],
                                        Uk_dense_row->x[ci]));
                SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[pi], U->v[i]->x[pi],
                                        Uk_dense_row->x[Q[k]]));
                num_of_fillin --;
                h[ci] = -1;
            }
        }

        // CPPU shortcuts
        pi = Ucx[pks];  // pointer for U(i,Q(ks))
        ASSERT(U->v[i]->i[pi] == Q[ks]);
        int64_t Ui_nz = U->v[i]->nz;
        if (i != ks)
        {
            // Mathematically, we should flip sign for U(i, Q(ks)). However,
            // since column ks of U will become column k after permutation,
            // which will be deleted when finished, we will instead delete
            // U(i,Q(ks)).
            if (num_of_fillin > 0)
            {
                // store U(i,Q(ks)) if there is fillin to be added
                SPEX_CHECK(SPEX_mpz_neg(Uiks, U->v[i]->x[pi]));
            }
            Ui_nz--;
            U->v[i]->i[pi] = U->v[i]->i[Ui_nz];
            SPEX_CHECK(SPEX_mpz_swap(U->v[i]->x[pi], U->v[i]->x[Ui_nz]));
            // skip scaling for L(:,i) by setting S(1,i) = S(1,i)*pending_scale
            SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 1, i),
                                    SPEX_LUU_2D(S, 1, i), pending_scale));
            // sd[i] = sd[i]*pending_scale
            SPEX_CHECK(SPEX_mpz_divexact(sd[i], sd[i],
                                    SPEX_MPQ_DEN(pending_scale)));
            SPEX_CHECK(SPEX_mpz_mul(sd[i], sd[i],
                                    SPEX_MPQ_NUM(pending_scale)));
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_neg(U->v[i]->x[pi], U->v[i]->x[pi]));
            SPEX_CHECK(SPEX_mpq_neg(SPEX_LUU_2D(S, 1, i), SPEX_LUU_2D(S, 1, i)));
            SPEX_CHECK(SPEX_mpz_neg(sd[i], sd[i]));
        }

        // finish RwSOP by checking if there is fillin that should be added
        if (num_of_fillin > 0)
        {
            // allocate additional space if needed
            if (U->v[i]->nzmax < Ui_nz+num_of_fillin)
            {
                SPEX_CHECK(SPEX_vector_realloc(U->v[i], Ui_nz+num_of_fillin));
            }
            // add FILLIN
            for (pk = 0; pk < Uk_dense_row->nz; pk++)
            {
                ck = Uk_dense_row->i[pk];
                if (h[ck] == -2)
                {
                    U->v[i]->i[Ui_nz] = ck;
                    // U(i,ck)= -U(i,Q(ks))*U(k,ck)/U(k,Q(k))
                    //        = Uiks       *U(k,ck)/U(k,Q(k))
                    if (i != ks)
                    {
                        SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[Ui_nz],
                                Uiks      , Uk_dense_row->x[ck]));
                    }
                    else
                    {
                        // sign has been flipped
                        SPEX_CHECK(SPEX_mpz_mul(U->v[i]->x[Ui_nz],
                            U->v[i]->x[pi], Uk_dense_row->x[ck]));
                    }
                    SPEX_CHECK(SPEX_mpz_divexact(U->v[i]->x[Ui_nz],
                        U->v[i]->x[Ui_nz], Uk_dense_row->x[Q[k]]));
                    Ui_nz++;
                    h[ck] = -1;
                    num_of_fillin --;
                    if (num_of_fillin == 0) {break;}
                }
            }
        }
        U->v[i]->nz = Ui_nz;
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
    Uk_nz = 0;
    for (pk = 0; pk < Uk_dense_row->nz; pk++)
    {
        j = Uk_dense_row->i[pk];
        if (j == Q[k])
        {
            // no need to copy U(k,Q[k]) since column k will be deleted anyway
            SPEX_CHECK(SPEX_mpz_set_ui(Uk_dense_row->x[j], 0));
            continue;
        }

        SPEX_CHECK(SPEX_mpz_swap(U->v[k]->x[Uk_nz], Uk_dense_row->x[j]));
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
            SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[ks]->x[pks]));
            if (sgn == 0) {continue;}

            if (Q_inv[j] == ks)
            {
                // this entry should be considered as the IPGE update of the
                // entry in the Q[k]-th column
                j = Q[k];
            }
            else if (Q_inv[j] < *jnext)
            {
                // update the index for the next nnz in current row
                *jnext = Q_inv[j];
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
    int64_t tmpi;
    tmpi = Q[k];      Q[k] = Q[ks];          Q[ks] = tmpi;
    Q_inv[Q[k]] = k;  Q_inv[Q[ks]] = ks;
    GOTCHA;

    //-------------------------------------------------------------------------
    // flip sign for columns and rows ks+1 to n and update corresponding sd
    //-------------------------------------------------------------------------
    for (i = ks+1; i < n; i++)
    {
        SPEX_CHECK(SPEX_mpq_neg(SPEX_LUU_2D(S, 3, i), SPEX_LUU_2D(S, 3, i)));
        SPEX_CHECK(SPEX_mpz_neg(sd[i], sd[i]));
    }
    GOTCHA;

    SPEX_FREE_ALL;
    return SPEX_OK;
}
