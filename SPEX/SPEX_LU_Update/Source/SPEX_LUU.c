//------------------------------------------------------------------------------
// SPEX_CHOLMOD/SPEX_LUU.c: perform LU update for column replacement
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to perform LU update for column replacement.
// L and U are modified regardless of success or failure.

#define SPEX_FREE_ALL                \
    SPEX_FREE(h);                    \
    SPEX_FREE(h_for_vk);             \
    SPEX_FREE(map);                  \
    SPEX_FREE(Ldiag);                \
    SPEX_FREE(Lr_offdiag);           \
    SPEX_FREE(Uci);                  \
    SPEX_FREE(Ucp);                  \
    SPEX_FREE(Ucx);                  \
    spex_scattered_vector_free(&Lk_dense_col);\
    spex_scattered_vector_free(&Uk_dense_row);\
    spex_scattered_vector_free(&vk_dense);    \
    SPEX_MPQ_CLEAR(one);

#include "spex_lu_update_internal.h"

SPEX_info SPEX_LUU
(
    SPEX_mat *A,         // the original matrix in compressed-column form
    SPEX_mat *L,         // stored in compressed-column form
    SPEX_mat *U,         // stored in compressed-row form
    mpz_t *d,               // an array of size n that stores the unscaled pivot
    mpz_t *sd,              // an array of size n that stores the scaled pivot
    mpq_t *S,               // an array of size 3*n that stores pending scales
    int64_t *P,             // row permutation
    int64_t *P_inv,         // inverse of row permutation
    int64_t *Q,             // column permutation
    int64_t *Q_inv,         // inverse of column permutation
    SPEX_vector **vk,       // pointer to the inserted column, which will be
                            // swapped with A->v[k] in the output if succeed
    int64_t k,              // the column index that vk will be inserted
    const SPEX_options *option// command parameters
)
{
    // initialize workspace
    SPEX_info info;
    if (!spex_initialized()) {return SPEX_PANIC;}

    if (L->n != U->n || L->m != L->n || U->m != U->n)
    {
        return SPEX_INCORRECT_INPUT;
    }
    int sgn_vkk, sgn_vkn, r;
    int64_t ks, p, i, j, inext, jnext, n = L->n;
    int64_t *h = NULL, *h_for_vk = NULL, *Ldiag = NULL, *Lr_offdiag = NULL,
        *Uci = NULL, *Ucp = NULL, *Ucx = NULL, *map = NULL;
    spex_scattered_vector *Lk_dense_col = NULL, *Uk_dense_row = NULL,
        *vk_dense = NULL;
    mpq_t one; SPEX_MPQ_SET_NULL(one);
    SPEX_CHECK(SPEX_mpq_init(one));
    SPEX_CHECK(SPEX_mpq_set_ui(one, 1, 1));

    h        = (int64_t*) SPEX_calloc(n, sizeof(int64_t));
    h_for_vk = (int64_t*) SPEX_calloc(n, sizeof(int64_t));
    map      = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    if (!h || !h_for_vk || !map)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    // get the row-wise nnz pattern for L and column-wise nnz pattern for U
    SPEX_CHECK(spex_get_nnz_pattern(&Ldiag, &Lr_offdiag, &Uci, &Ucp, &Ucx,
        L, U, P, option));

    // initialize environment for the inserted column
    SPEX_CHECK(spex_get_scattered_v(&vk_dense, *vk, n, true));
    //index of column of L used to perform IPGE for vk
    int64_t last_update = -1;
    // the index of last nnz in vk(P[last_update+1,n-2])
    int64_t vk_2ndlastnz = -1;
    int64_t maximum, last_max_ks;
    int use_col_n = 0; // 0: unknown; 1: use; -1: don't use
    for (p = 0; p < vk_dense->nz; p++)
    {
        i = vk_dense->i[p];
        if (P_inv[i] > vk_2ndlastnz && P_inv[i] != n-1)
        {
            vk_2ndlastnz = P_inv[i];
        }
        h_for_vk[i] = -1;
    }

    // update column k of A by swapping A->v[k] with vk
    SPEX_vector *tmpv;
    tmpv = *vk; *vk = A->v[k]; A->v[k] = tmpv;

    k = Q_inv[k];
    last_max_ks = k;
    // remove entries in column k of U, but ignore the diagnal of k-th row of U
#ifdef SPEX_DEBUG
    printf("k=%ld;\nuk=[",k);
    j = 0;
#endif
    for (p = Ucp[Q[k]]; p < Ucp[Q[k]+1]-1; p++)
    {
        i = Uci[p];
        if (i < k)
        {
#ifdef SPEX_DEBUG
            while(j < i) { printf("0 ");j++;}
            SPEX_CHECK(SPEX_gmp_printf("%Zd ",U->v[i]->x[Ucx[p]]));
            j++;
#endif
            // move the last entry to current position
            U->v[i]->nz--;
            SPEX_CHECK(SPEX_mpz_swap(U->v[i]->x[Ucx[p]],
                                     U->v[i]->x[U->v[i]->nz]));
            U->v[i]->i[Ucx[p]] = U->v[i]->i[U->v[i]->nz];
        }
    }
#ifdef SPEX_DEBUG
    for (; j < n; j++) { printf("0 ");}
    printf("]';\n");
#endif

    if (k < n-1)
    {
        // build Lk_dense_col and Uk_dense_row, remove explicit 0
        SPEX_CHECK(spex_get_scattered_v(&Lk_dense_col, L->v[k], n, false));
        SPEX_CHECK(spex_get_scattered_v(&Uk_dense_row, U->v[k], n, false));

        // initialize certain variables required by the loop
        SPEX_CHECK(spex_find_next_nz(&inext, Lk_dense_col, P_inv, k));
        SPEX_CHECK(spex_find_next_nz(&jnext, Uk_dense_row, Q_inv, k));
    }

    // push column k to position n-1
    while (k < n-1)
    {
#if SPEX_DEBUG
        printf("-----------------------------------------------------------");
        printf("-------------\nl%ld=[\n",k);
        for (int64_t ii =0;ii<n;ii++)
        {
            if (ii==k)
            {
                for(int64_t pp = 0;pp<n;pp++)
                {
                    SPEX_CHECK(SPEX_gmp_printf("%Zd ",Lk_dense_col->x[pp])); 
                }
                for(int64_t pp=0;pp<L->v[k]->nzmax;pp++)
                {
                    SPEX_CHECK(SPEX_mpz_sgn(&r,L->v[k]->x[pp]));
                    ASSERT(r==0);
                }
            }
            else
            {
                for(int64_t jj=0;jj<n;jj++)
                {
                    bool found_the_entry = false;
                    for (int64_t pp = 0;pp<L->v[ii]->nz;pp++)
                    {
                        if (jj==L->v[ii]->i[pp])
                        {
                        SPEX_CHECK(SPEX_gmp_printf("%Zd ",L->v[ii]->x[pp]));
                        found_the_entry = true;
                        break;
                        }
                    }
                    if (!found_the_entry){printf("0 ");}
                }
            }
            printf("%s;\n",ii==n-1?"]":"");
        }
        printf("\nl%ld=l%ld'*diag([",k,k);
        for(int64_t ii = 0; ii<n;ii++)
        {
            SPEX_CHECK(SPEX_gmp_printf("%Qd*%Qd ", SPEX_LUU_2D(S, 1, ii),
                                                   SPEX_LUU_2D(S,3,ii)));
        }
        printf("]);\nu%ld=[\n",k);
        for (int64_t ii =0;ii<n;ii++)
        {
            if (ii==k)
            {
                for(int64_t pp = 0;pp<n;pp++)
                {
                    SPEX_CHECK(SPEX_gmp_printf("%Zd ",Uk_dense_row->x[pp])); 
                }
                for(int64_t pp=0;pp<U->v[k]->nzmax;pp++)
                {
                    SPEX_CHECK(SPEX_mpz_sgn(&r,U->v[k]->x[pp]));
                    ASSERT(r==0);
                }
            }
            else
            {
                for(int64_t jj=0;jj<n;jj++)
                {
                    bool found_the_entry = false;
                    for (int64_t pp = 0;pp<U->v[ii]->nz;pp++)
                    {
                        if (jj==U->v[ii]->i[pp])
                        {
                        SPEX_CHECK(SPEX_gmp_printf("%Zd ",U->v[ii]->x[pp]));
                        found_the_entry = true;
                        break;
                        }
                    }
                    if (!found_the_entry){printf("0 ");}
                }
            }
            printf("%s;\n",ii==n-1?"]":"");
        }
        printf("\nuk = u%ld(:,k)+uk;\n u%ld(:,k) = uk;\nu%ld=diag([",k,k,k);
        for(int64_t ii = 0; ii<n;ii++)
        {
            SPEX_CHECK(SPEX_gmp_printf("%Qd*%Qd ", SPEX_LUU_2D(S, 2, ii),
                                                   SPEX_LUU_2D(S, 3, ii)));
        }
        printf("])*u%ld;\nsd%ld=[",k,k);
        for(int64_t ii=0;ii<n;ii++) SPEX_CHECK(SPEX_gmp_printf("%Zd ",sd[ii]));
        printf("];\nd%ld=[",k);
        for(int64_t ii=0;ii<n;ii++) SPEX_CHECK(SPEX_gmp_printf("%Zd ",d[ii]));
        printf("];\nD%ld=diag(1./(sd%ld'.*[1;sd%ld(1:end-1)']));\n",k,k,k);
        printf("difference%ld = norm(A0-l%ld*D%ld*u%ld)\n",k,k,k,k);

            printf("\nh=[");
        for(int64_t ii=0;ii<n;ii++)printf("%ld ",h[ii]);
            printf("]\nh4vk=[");
        for(int64_t ii=0;ii<n;ii++)printf("%ld ",h_for_vk[ii]);
            printf("]\nvk_nnz pattern(%ld)=[",vk_dense->nz);
        for(int64_t ii=0;ii<vk_dense->nz;ii++)printf("%ld ",vk_dense->i[ii]);
            printf("]\nP=[");
        for(int64_t ii=0;ii<n;ii++)printf("%ld ",P[ii]);
        printf("]\nQ=[");
        for(int64_t ii=0;ii<n;ii++)printf("%ld ",Q[ii]);
        printf("]\nuse_col_n=%d\n",use_col_n);
#endif

        // no need to update vk if we know not to use it
        if (use_col_n >= 0)
        {
            // get the (k-1)-th IPGE update of inserted column, so that vk
            // will be ready to be inserted as column k if needed
            SPEX_CHECK(spex_triangular_solve(vk_dense, h_for_vk,
                &last_update, &vk_2ndlastnz, k, L, U, Ldiag, Ucp, Ucx,
                S, (const mpz_t*)sd, d, P, P_inv, Q));
#ifdef SPEX_DEBUG
            printf("vk=[ ");
            for (int64_t ii=0;ii<vk_dense->nzmax;ii++)
            {
                SPEX_CHECK(SPEX_gmp_printf("%Zd  ",vk_dense->x[ii]));
            }
            printf("];\n");
#endif
        }
        SPEX_CHECK(SPEX_mpz_sgn(&sgn_vkk, vk_dense->x[P[k]]));
        SPEX_CHECK(SPEX_mpz_sgn(&sgn_vkn, vk_dense->x[P[n-1]]));

        if (jnext < n)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&r, Uk_dense_row->x[Q[jnext]]));
            if (r == 0)
            {
                SPEX_CHECK(spex_find_next_nz(&jnext, Uk_dense_row, Q_inv, k));
            }
        }

        // report singular if 
        // - remaining entries in current row of U are 0s and the current row
        //   of vk is also 0.
        // - OR all entries below (k-1)-th row in vk are zeros
        // - OR all off-diagonal entries below (k-1)-th row in vk
        //   and column n-1 of U are zeros
        // That is, report singular for any of the following case
        // x 0 0 0 0 0          x . . . . 0          x . . . 0 0   <- row k
        // . x . . . .          . x . . . 0          . x . . 0 0
        // . . x . . .    or    . . x . . 0    or    . . x . 0 0
        // . . . x . .          . . . x . 0          . . . x 0 0
        // . . . . x .          . . . . x 0          . . . . x x
        //           ^                    ^                    ^
        //           |                    |                    |
        //          vk                    vk                   vk
        if ((jnext == n && sgn_vkk == 0 ) ||
            (vk_2ndlastnz == -1 &&
             (sgn_vkn == 0 ||
              Ucp[Q[n-1]+1]-Ucp[Q[n-1]] == 1 || Uci[Ucp[Q[n-1]+1]-2] < k)))
        {
            SPEX_FREE_ALL;
            return SPEX_SINGULAR;
        }
        // if the next nnz in current row is in vk, then use vk, which will
        // help to avoid performing extra ipge iterations for vk
        if (jnext == n && use_col_n >= 0)
        {
            // if jnext == n, swapping columns k and vk will be more efficient,
            // since there is no need to backtrack vk (i.e., col n) (we just
            // perform k-th IPGE iteration for column n), and column n-1 can be
            // updated by scaling after CPPU.
            ks = n;
            SPEX_CHECK(spex_finalize_and_insert_vk(vk_dense, h_for_vk,
                U, L, S, d, Ldiag, (const mpz_t*)sd, Q, P_inv, k, k, one));

            SPEX_CHECK(spex_cppu(L, U, S, d, sd, Lk_dense_col,
                Uk_dense_row, &inext, &jnext, h, Q, Q_inv,
                P, P_inv, Ldiag, Uci, Ucp, Ucx, k, ks));
            break;
        }
        if (inext < n)
        {
            ASSERT (inext > k);
            SPEX_CHECK(SPEX_mpz_sgn(&r, Lk_dense_col->x[P[inext]]));
            if (r == 0)
            {
                SPEX_CHECK(spex_find_next_nz(&inext, Lk_dense_col, P_inv, k));
            }
        }
        printf("k(%ld) inext(%ld) jnext(%ld)\n",k, inext,jnext);

        //----------------------------------------------------------------------
        // if L(:,k) has zero off-diagonal, then only perform dppu, which will
        // maintain the sparsity of L(:,k). Use dppu1 if possible.
        // When arriving the last iteration, always use the inserted column
        // if possible, since we can perform less IPGE iterations for it.
        //----------------------------------------------------------------------
        if (inext == n)
        {
            // force S(3,k) = S(2,k)*S(3,k) and S(2,k) = 1 since we only care
            // about the row in frame k and simply treat the column as 0
            SPEX_CHECK(SPEX_mpq_equal(&r, SPEX_LUU_2D(S, 2, k), one));
            if (r == 0) //S(2,k) != 1
            {
                SPEX_CHECK(SPEX_mpq_mul(SPEX_LUU_2D(S, 3, k),
                                        SPEX_LUU_2D(S, 3, k), SPEX_LUU_2D(S, 2, k)));
                SPEX_CHECK(SPEX_mpq_set(SPEX_LUU_2D(S, 2, k), one));
            }

            // build the map to find ks if current map is out of date.
            // all the swaps (i.e., pivot updates) except the last one
            // using this map will not change the nnz patter of current
            // frame, since only scaling will be involved.
            // x 0 0 0 0 0 0 0 0 x 
            // 0 x . 0 . . . . . .
            // 0 . x 0 . . . . . .
            // 0 0 0 x 0 . . . . .
            // 0 . . . x 0 0 0 0 .
            // 0 . . . 0 x . . 0 .
            // 0 . . . 0 . x . 0 .
            // 0 . . . 0 . . x 0 .
            // 0 . . . . 0 0 0 x .
            // 0 . . . . . . . . x
            // ^     ^ ^       ^
            // |     | |       |
            // 2     5 6       10
            // we will swap 2 and 5 first, then with 6, and with 10
            // sequencially, which can be found from the map eazily
            if (n-2 > last_max_ks)
            {
                for (j = k+1; j <= n-2; j++)
                {
                    if (Ucp[Q[j]+1] - Ucp[Q[j]] == 1)
                    {
                        // no off-diag in column Q[j] of U
                        maximum = SPEX_MAX(Lr_offdiag[P[j]], -1);
                    }
                    else
                    {
                        maximum = SPEX_MAX(Lr_offdiag[P[j]],Uci[Ucp[Q[j]+1]-2]);
                    }
                    for (i = k; i < j;)
                    {
                        if (maximum <= i)
                        {
                            map[i] = j;
                            break;
                        }
                        i = map[i];
                    }
                }
                last_max_ks = n-2;
            }

            // use the inserted column only when its last entry is nnz and
            // using it instead of column n-1 can make a bigger jump.
            if (use_col_n == 0 && Lr_offdiag[P[n-1]] <= k)
            {
                if (sgn_vkn == 0)// vk_dense[P[n-1]] == 0
                {
                    // the inserted column cannot be used
                    use_col_n = -1;
                }
                else
                {
                    // use vk only when the index of off diagonal entry in
                    // column Q(n-1) of U is larger than vk_2ndlastnz, that is
                    //        x . . . . .  <- row k
                    //        0 x . . . 0
                    //        0 . x . x 0
                    //        0 . . x 0 0
                    //        0 0 0 0 x x
                    //                  ^
                    //                  |
                    //                 vk at (k-1)-th IPGE
                    if (Ucp[Q[n-1]+1] - Ucp[Q[n-1]] > 1 && vk_2ndlastnz <= k &&
                        Uci[Ucp[Q[n-1]+1]-2] > vk_2ndlastnz)
                    {
                        use_col_n = 1;
                    }
                    else
                    {
                        use_col_n = -1;
                    }
                }
                if (use_col_n == -1)
                {
                    if (Ucp[Q[n-1]+1] - Ucp[Q[n-1]] == 1)
                    {
                        // no off-diag in column Q[n-1] of U
                        maximum = SPEX_MAX(Lr_offdiag[P[n-1]], -1);
                    }
                    else
                    {
                        maximum = SPEX_MAX(Lr_offdiag[P[n-1]],
                                           Uci[Ucp[Q[n-1]+1]-2]);
                    }
                }
                else
                {
                    maximum = SPEX_MAX(Lr_offdiag[P[n-1]], vk_2ndlastnz);
                }
                for (i = k; i < n-1;)
                {
                    if (maximum <= i)
                    {
                        map[i] = n-1;
                        break;
                    }
                    i = map[i];
                }

                last_max_ks = n-1;
            }

            // get ks from the map
            ks = map[k];
            ASSERT(ks > k);
            if (ks == n-1 && use_col_n == 1)
            {
                SPEX_CHECK(spex_finalize_and_insert_vk(vk_dense, h_for_vk, U,
                    L, S, d, Ldiag, (const mpz_t*)sd, Q, P_inv, k, n-1, one));
                ks = n;
            }
            if (jnext > ks || (ks == n && jnext >= n-1 && sgn_vkk == 0))
            {
                SPEX_CHECK(spex_dppu1(L, U, S, d, sd, Lk_dense_col,
                    Uk_dense_row, &inext, h, Q, Q_inv, P, P_inv,
                    Ldiag, Uci, Ucp, Ucx, k, ks));
            }
            else
            {
                SPEX_CHECK(spex_dppu2(L, U, S, d, sd, Lk_dense_col,
                    Uk_dense_row, &jnext, h, Q, Q_inv, P, P_inv,
                    Ldiag, Uci, Ucp, Ucx, k, ks));
            }
        }
        else
        {
            // The case when jnext == n (that is, the remaining entries in the
            // current row except the n-th column (vk) are all zeroes) is
            // handled. And the case when all remaining entries in current
            // row are zeroes will cause singularity. Therefore, the cases
            // left to consider here are U(k,Q(n-1:n))=[x,0] and
            // U(k,Q(n-1:n)) = [x,x]
            if (jnext == n-1)
            {
                if (sgn_vkk == 0) // U(k,n) = vk[P[k]] == 0
                {
                    // if jnext == n-1 and U(k,n) == 0, use dppu1 with column n
                    // only when we see the following pattern
                    //
                    //             x 0 0 0 . 0
                    //             0 x . . . 0
                    //             0 . x . . 0
                    //             0 . . x . 0
                    //             . 0 0 0 x x
                    //                       ^
                    //                       |
                    //                   col n or vk
                    //
                    if (vk_2ndlastnz == -1/*then vk[P[n-1]] must != 0*/  &&
                        Lr_offdiag[P[n-1]] <= k && inext >= n-1)
                    {
                        // perform diagnal swapping with columns n
                        ks = n;
                        SPEX_CHECK(spex_finalize_and_insert_vk(vk_dense,
                            h_for_vk, U, L, S, d, Ldiag, (const mpz_t*)sd, Q,
                            P_inv, k, n-1, one));
                        SPEX_CHECK(spex_dppu1(L, U, S, d, sd, Lk_dense_col,
                            Uk_dense_row, &inext, h, Q, Q_inv, P,
                            P_inv, Ldiag, Uci, Ucp, Ucx, k, ks));
                        break;
                    }
                }
                else
                {
                    // if jnext == n-1 and U(k,n) != 0, swapping columns k and
                    // n will be more efficient, since no matter column n or
                    // n-1 is swapped, column n needs n IPGE iterations, while
                    // if column n-1 is swapped, additional backtracking for
                    // column n-1 needs to be performed.
                    ks = n;
                    SPEX_CHECK(spex_finalize_and_insert_vk(vk_dense,
                        h_for_vk, U, L, S, d, Ldiag, (const mpz_t*)sd, Q,
                        P_inv, k, k, one));
                    SPEX_CHECK(spex_cppu(L, U, S, d, sd, Lk_dense_col,
                        Uk_dense_row, &inext, &jnext, h, Q, Q_inv, P, P_inv,
                        Ldiag, Uci, Ucp, Ucx, k, ks));
                    break;
                }
            }
            if (inext == k+1 ||
                (Ucp[Q[jnext]+1] - Ucp[Q[jnext]] >= 2 &&
                 Uci[Ucp[Q[jnext]+1]-2] <= k /*Since Uci,Ucp are nnz pattern of
                 U found ahead, and potential fillins added to row k of U
                 during the update process are not included, the 3rd condition
                 is "<=" instead of "==". Uk_dense_row->x[Q[jnext]] can be
                 double checked to make sure that it is nnz*/))
            {
                // use cppu if we see one of the following patterns
                // x . . . .            x 0 0 0 x
                // x x . . .            . x . . 0
                // . . x . .     or     . . x . 0
                // . . . x .            . . . x 0
                // . . . . x            . . . . x
                // These implicitly include the case of jnext == k+1.
                // jnext < n holds since jnext == n has been handled,
                ks = jnext;
                SPEX_CHECK(spex_cppu(L, U, S, d, sd, Lk_dense_col, Uk_dense_row,
                    &inext, &jnext, h, Q, Q_inv, P, P_inv, Ldiag,
                    Uci, Ucp, Ucx, k, ks));
            }
            else
            {
                ks = SPEX_MIN(n-2, (inext < jnext) ? inext: jnext-1);
                // build the map to find ks if current map is out of date.
                // all the swaps (i.e., pivot updates) except the last one
                // using this map will not change the nnz patter of current
                // frame, since only scaling will be involved.
                // x 0 0 0 0 0 0 0 0 x 
                // 0 x . 0 . . . . . .
                // 0 . x 0 . . . . . .
                // 0 0 0 x 0 . . . . .
                // 0 . . . x 0 0 0 0 .
                // 0 . . . 0 x . . 0 .
                // 0 . . . 0 . x . 0 .
                // 0 . . . 0 . . x 0 .
                // 0 . . . . 0 0 0 x .
                // 0 . . . . . . . . x
                // ^     ^ ^       ^
                // |     | |       |
                // 2     5 6       10
                // we will swap 2 and 5 first, then with 6, and with 10
                // sequencially, which can be found from the map eazily
                if (ks > last_max_ks)
                {
                    ASSERT(k+1 >= last_max_ks);
                    for (j = k+1; j <= ks; j++)
                    {
                        if (Ucp[Q[j]+1] - Ucp[Q[j]] == 1)
                        {
                            // no off-diag in column Q[j] of U
                            maximum = SPEX_MAX(Lr_offdiag[P[j]], -1);
                        }
                        else
                        {
                            maximum = SPEX_MAX(Lr_offdiag[P[j]],
                                               Uci[Ucp[Q[j]+1]-2]);
                        }
                        for (i = k; i < j;)
                        {
                            if (maximum <= i)
                            {
                                map[i] = j;
                                break;
                            }
                            i = map[i];
                        }
                    }
                    last_max_ks = ks;
                }

                // get ks from the map
                ks = map[k];

                SPEX_CHECK(spex_dppu1(L, U, S, d, sd, Lk_dense_col,
                    Uk_dense_row, &inext, h, Q, Q_inv, P, P_inv,
                    Ldiag, Uci, Ucp, Ucx, k, ks));
            }
        }

        // update k
        ASSERT(ks > k);
        if(ks != n)
        {
            k = ks;
        }
        else
        {
            break;
        }
    }

    // k will be the column index where vk is inserted
    if (k == n-1)
    {
        SPEX_CHECK(spex_triangular_solve(vk_dense, h_for_vk, 
            &last_update, NULL /*&vk_2ndlastnz*/, k, L, U, Ldiag, Ucp, Ucx,
            S, (const mpz_t*)sd, d, P, P_inv, Q));
#ifdef SPEX_DEBUG
            printf("vk=[ ");
            for (int64_t ii=0;ii<vk_dense->nzmax;ii++)
            {
                SPEX_CHECK(SPEX_gmp_printf("%Zd  ",vk_dense->x[ii]));
            }
            printf("];\n");
#endif
        // check again in case k is initially n-1
        SPEX_CHECK(SPEX_mpz_sgn(&sgn_vkn, vk_dense->x[P[n-1]]));
        if (sgn_vkn == 0)
        {
            SPEX_FREE_ALL;
            return SPEX_SINGULAR;
        }
        SPEX_CHECK(spex_finalize_and_insert_vk(vk_dense, h_for_vk, U, L, S, d,
            Ldiag, (const mpz_t*)sd, Q, P_inv, k, k, one));
        // U(n-1,n-1) = d[n-1]
        SPEX_CHECK(SPEX_mpz_set(U->v[n-1]->x[0], d[n-1]));
        U->v[n-1]->i[0] = Q[n-1];
        U->v[n-1]->nz = 1;
        // update d[k]=vk[k], sd[k]=U(k,k)=vk[k]
        SPEX_CHECK(SPEX_mpz_set(sd[n-1], d[n-1]));
    }
    else
    {
        // update U(k,Q(k)) and S(:,k)
        if (U->v[k]->nzmax <= U->v[k]->nz)
        {
            SPEX_CHECK(SPEX_vector_realloc(U->v[k], U->v[k]->nzmax+1));
        }
        SPEX_CHECK(SPEX_mpz_set(U->v[k]->x[U->v[k]->nz], sd[k]));
        U->v[k]->i[U->v[k]->nz] = Q[k];
        U->v[k]->nz++;
    }
    // S(:,k)=[1;1;1]
    SPEX_CHECK(SPEX_mpq_set_ui(SPEX_LUU_2D(S, 1, k), 1, 1));
    SPEX_CHECK(SPEX_mpq_set_ui(SPEX_LUU_2D(S, 2, k), 1, 1));
    SPEX_CHECK(SPEX_mpq_set_ui(SPEX_LUU_2D(S, 3, k), 1, 1));

//#ifdef SPEX_DEBUG
    // check if A=LD^(-1)U
    bool result;
    SPEX_CHECK(spex_verify(&result, L, U, A, h, (const mpz_t*) sd, d,
        S, P, P_inv, Q, Q_inv, Ldiag, Ucp, Ucx, option));
    printf("the factorization is %s\n", result?"correct":"incorrect");
//#endif

    SPEX_FREE_ALL;
    return SPEX_OK;
    //return result?SPEX_OK:SPEX_INCORRECT_INPUT;
}
