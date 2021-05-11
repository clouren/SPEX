//------------------------------------------------------------------------------
// SPEX_Update/SPEX_Update_LU_ColRep.c: perform LU update for column replacement
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function is to perform LU update for column replacement.
// L and U are modified regardless of success or failure.
// In the output, vk will be swapped with A->v[k].

#define SPEX_FREE_ALL                \
    SPEX_FREE(h);                    \
    SPEX_FREE(h_for_vk);             \
    SPEX_FREE(map);                  \
    SPEX_FREE(Lr_offdiag);           \
    SPEX_FREE(Uc_offdiag);           \
    SPEX_FREE(Uci);                  \
    SPEX_FREE(Ucx);                  \
    spex_scattered_vector_free(&Lk_dense_col);\
    spex_scattered_vector_free(&Uk_dense_row);\
    spex_scattered_vector_free(&vk_dense);    \
    SPEX_MPQ_CLEAR(one);

#include "spex_update_internal.h"

#define SL(k) (L->v[(k)]->scale)
#define SU(k) (U->v[(k)]->scale)

SPEX_info SPEX_Update_LU_ColRep//SPEX_colrep_ modcol
(
    SPEX_mat *A,            // the original matrix in compressed-column form
    SPEX_mat *L,            // stored in compressed-column form
    SPEX_mat *U,            // stored in compressed-row form
    mpz_t *sd,              // an array of size n that stores the scaled pivot
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
    int64_t *h = NULL, *h_for_vk = NULL, *Lr_offdiag = NULL, *Uc_offdiag = NULL,
        *Uci = NULL, *Ucx = NULL, *map = NULL;
    spex_scattered_vector *Lk_dense_col = NULL, *Uk_dense_row = NULL,
        *vk_dense = NULL;
    mpq_t one; SPEX_MPQ_SET_NULL(one);
    SPEX_CHECK(SPEX_mpq_init(one));
    SPEX_CHECK(SPEX_mpq_set_ui(one, 1, 1));

    h          = (int64_t*) SPEX_calloc(n, sizeof(int64_t));
    h_for_vk   = (int64_t*) SPEX_calloc(n, sizeof(int64_t));
    Uc_offdiag = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    Lr_offdiag = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    Uci        = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    Ucx        = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    map        = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    if (!h || !h_for_vk || !Uci || !Ucx || !map)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    // initialize environment for the inserted column
    SPEX_CHECK(spex_update_get_scattered_v(&vk_dense, *vk, n, true));
    //index of column of L used to perform IPGE for vk
    int64_t last_update = -1;
    // the index of last nnz in vk(P[last_update+1,n-2])
    int64_t vk_2ndlastnz = -1;
    int64_t maximum, last_max_ks, LU_scaned, Uc_jnext_nz = 0;
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

    // initialize for the while loop
    LU_scaned = Q_inv[k];
    // remove column k of U
    for (i = 0; i < LU_scaned; i++)
    {
        for (p = 1 ; p < U->v[i]->nz; p++)
        {
            j = U->v[i]->i[p];
            if (j == k)
            {
                // move the last entry to current position
                U->v[i]->nz--;
                SPEX_CHECK(SPEX_mpz_swap(U->v[i]->x[p],
                                         U->v[i]->x[U->v[i]->nz]));
                U->v[i]->i[p] = U->v[i]->i[U->v[i]->nz];
                break;
            }
        }
    }
    k = LU_scaned; // = Q_inv[k];
    last_max_ks = k;
    // initialize entries k+1:n
    for (i = k+1; i < n ; i++)
    {
        Lr_offdiag[i] = -1;
        Uc_offdiag[i] = -1;
    }

    if (k < n-1)
    {
        // build Lk_dense_col and Uk_dense_row, remove explicit 0
        SPEX_CHECK(spex_update_get_scattered_v(&Lk_dense_col, L->v[k], n, false));
        SPEX_CHECK(spex_update_get_scattered_v(&Uk_dense_row, U->v[k], n, false));

        // initialize certain variables required by the loop
        SPEX_CHECK(spex_update_find_next_nz(&inext, Lk_dense_col, P_inv, k));
        SPEX_CHECK(spex_update_find_next_nz(&jnext, Uk_dense_row, Q_inv, k));
    }

    // push column k to position n-1
    while (k < n-1)
    {
#ifdef SPEX_DEBUG
        //printf("-----------------------------------------------------------");
        //printf("-------------\n\n");
        for (int64_t ii = 0; ii < n; ii++)
        {
            if (ii != k && L->v[ii]->i[0] != P[ii])
            {
                printf("Incorrect col %ld of L\n",ii);
                SPEX_CHECK(SPEX_PANIC);
            }
            SPEX_CHECK(SPEX_mpz_sgn(&r, Lk_dense_col->x[ii]));
            if (r != 0)
            {
                bool gotcha = false;
                for (int64_t pp = 0; pp < Lk_dense_col->nz; pp++)
                {
                    if (Lk_dense_col->i[pp] == ii)
                    {
                        gotcha = true;
                        break;
                    }
                }
                if (gotcha == false)
                {
                    printf("file %s line %d\n",__FILE__,__LINE__);
                    for (int64_t pp = 0; pp < Lk_dense_col->nz; pp++)
                    {
                        printf("%ld ", Lk_dense_col->i[pp]);
                    }
                    SPEX_CHECK(SPEX_gmp_printf("\nLk(%ld)=%Zd\n",
                        ii, Lk_dense_col->x[ii]));
                    SPEX_CHECK(SPEX_PANIC);
                }
            }
        }
        for (int64_t ii =0;ii<n;ii++)
        {
            if (ii != k && U->v[ii]->i[0] != Q[ii])
            {
                printf("Incorrect row %ld of U\n",ii);
                SPEX_CHECK(SPEX_PANIC);
            }
            SPEX_CHECK(SPEX_mpz_sgn(&r, Uk_dense_row->x[ii]));
            if (r != 0)
            {
                bool gotcha = false;
                for(int64_t pp = 0; pp < Uk_dense_row->nz; pp++)
                {
                    if (Uk_dense_row->i[pp] == ii)
                    {
                        gotcha = true;
                        break;
                    }
                }
                if (gotcha == false)
                {
                    printf("file %s line %d\n",__FILE__,__LINE__);
                    for (int64_t pp = 0; pp < Uk_dense_row->nz; pp++)
                    {
                        printf("%ld ", Uk_dense_row->i[pp]);
                    }
                    SPEX_CHECK(SPEX_gmp_printf("\nUk(%ld)=%Zd\n",
                        ii, Uk_dense_row->x[ii]));
                    SPEX_CHECK(SPEX_PANIC);
                }
            }
        }
        for (int64_t ii = 0; ii < n; ii++)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&r, sd[ii]));
            if (r == 0)
            {
                printf("sd[%ld]=0 at file %s line %d\n",ii,__FILE__,__LINE__);
                SPEX_CHECK(SPEX_PANIC);
            }
            if (h[ii] < -1)
            {
                printf("h[%ld]=%ld at file %s line %d\n", ii, h[ii],
                    __FILE__, __LINE__);
                SPEX_CHECK(SPEX_PANIC);
            }
        }
#endif

        // update vk if we might use it
        if (use_col_n >= 0)
        {
            // get the (k-1)-th IPGE update of inserted column, so that vk
            // will be ready to be inserted as column k if needed
            SPEX_CHECK(spex_update_triangular_solve(vk_dense, h_for_vk,
                &last_update, &vk_2ndlastnz, k, L, U,
                (const mpz_t*)sd, P, P_inv));
        }
        SPEX_CHECK(SPEX_mpz_sgn(&sgn_vkk, vk_dense->x[P[k]]));
        SPEX_CHECK(SPEX_mpz_sgn(&sgn_vkn, vk_dense->x[P[n-1]]));

        if (jnext < n)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&r, Uk_dense_row->x[Q[jnext]]));
            if (r == 0)
            {
                SPEX_CHECK(spex_update_find_next_nz(&jnext, Uk_dense_row, Q_inv, k));
            }
        }
        else // jnext == n, i.e., only 1 nnz in row k of U
        {
            ASSERT(Uk_dense_row->nz == 1);
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
        bool Is_Singular = false;
        // case 1
        if (jnext == n && sgn_vkk == 0 )
        {
            Is_Singular = true;
        }
        // case 2
        if (vk_2ndlastnz == -1 && sgn_vkn == 0)
        {
            Is_Singular = true;
        }
        // case 3, check only when whole U is scaned 
        if (vk_2ndlastnz == -1 && LU_scaned == n-1 && Uc_offdiag[n-1] <= k)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&r, Uk_dense_row->x[Q[n-1]]));
            if (r == 0)
            {
                Is_Singular = true;
            }
        }
        if (Is_Singular)
        {
            SPEX_gmp_printf("U(k,Q[n-1])=%Zd\n",Uk_dense_row->x[Q[n-1]]);
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
            SPEX_CHECK(spex_update_finalize_and_insert_vk(vk_dense, h_for_vk,
                U, L, (const mpz_t*)sd, Q, P_inv, k, k, one));

            SPEX_CHECK(spex_update_cppu(L, U, sd, Lk_dense_col,
                Uk_dense_row, &inext, &jnext, h, Q, Q_inv,
                P, P_inv, NULL, NULL, 0, k, ks));
            break;
        }
        if (inext < n)
        {
            ASSERT (inext > k);
            SPEX_CHECK(SPEX_mpz_sgn(&r, Lk_dense_col->x[P[inext]]));
            if (r == 0)
            {
                SPEX_CHECK(spex_update_find_next_nz(&inext, Lk_dense_col, P_inv, k));
            }
        }
#ifdef SPEX_DEBUG
        printf("k(%ld) inext(%ld) jnext(%ld)\n",k, inext,jnext);
#endif

        //----------------------------------------------------------------------
        // if L(:,k) has zero off-diagonal, then only perform dppu, which will
        // maintain the sparsity of L(:,k). Use dppu1 if possible.
        // When arriving the last iteration, always use the inserted column
        // if possible, since we can perform less IPGE iterations for it.
        //----------------------------------------------------------------------
        if (inext == n)
        {
            if (Lk_dense_col->nz != 1)
            {
                Lk_dense_col->nz = 1;
                Lk_dense_col->i[0] = P[k];
            }

            // scan L(:, k+1:n-1) and U(k+1:n-1, :)
            if (LU_scaned < n-1)
            {
                for (j = k+1; j < n; j++)
                {
#ifdef SPEX_DEBUG
                    // the first entry should be the corresponding pivot
                    ASSERT(P_inv[L->v[j]->i[0]] == j);
                    int sgn;
                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[j]->x[0]));
                    ASSERT (sgn != 0);
#endif
                    for (p = 1; p < L->v[j]->nz; p++)
                    {
                        // row index
                        i = P_inv[L->v[j]->i[p]];

                        // get the last row-wise off-diagonal entries
                        Lr_offdiag[i] = j; // get the column index
                    }
                }

                for (i = k+1; i < n; i++)
                {
                    ASSERT(Q_inv[U->v[i]->i[0]] == i);
                    for (p = 1 ; p < U->v[i]->nz; p++)
                    {
                        j = Q_inv[U->v[i]->i[p]];
                        Uc_offdiag[j] = i;
                    }
                }

                LU_scaned = n-1;
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
                    maximum = SPEX_MAX(Lr_offdiag[j], Uc_offdiag[j]);
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

            // Check if we should use the inserted column when either column
            // Q[n-1] or vk could be used
            //        x . . . . .          x . . . . .  <- row k
            //        0 x . . . 0          0 x . . 0 .
            //        0 . x . . 0    or    0 . x . 0 .
            //        0 . . x . 0          0 . . x 0 .
            //        0 0 0 0 x x          0 0 0 0 x .
            if (use_col_n == 0 && Lr_offdiag[n-1] <= k &&
                (vk_2ndlastnz <= k /*vk could be used*/ ||
                 Uc_offdiag[n-1] <= k /*col Q[n-1] could be used*/))
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
                    if (Uc_offdiag[n-1] > vk_2ndlastnz && vk_2ndlastnz <= k)
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
                    maximum = SPEX_MAX(Lr_offdiag[n-1], Uc_offdiag[n-1]);
                }
                else
                {
                    maximum = SPEX_MAX(Lr_offdiag[n-1], vk_2ndlastnz);
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
                SPEX_CHECK(spex_update_finalize_and_insert_vk(vk_dense,
                    h_for_vk, U, L, (const mpz_t*)sd, Q, P_inv, k, n-1, one));
                ks = n;
            }
            if (jnext > ks || (ks == n && jnext >= n-1 && sgn_vkk == 0))
            {
                SPEX_CHECK(spex_update_dppu1(L, U, sd, Lk_dense_col,
                    Uk_dense_row, &inext, h, Q, Q_inv, P, P_inv, k, ks));
            }
            else
            {
                SPEX_CHECK(spex_update_dppu2(L, U, sd, Lk_dense_col,
                    Uk_dense_row, &jnext, h, Q, Q_inv, P, P_inv, k, ks));
            }
        }
        else // if L(:,k) has more than 1 nnz
        {
            // scan L(:, k+1:jnext) and U(k+1:jnext, :)
            if (LU_scaned < jnext && (jnext > k+1 || jnext == n-1))
            {
                for (j = k+1; j <= jnext; j++)
                {
#ifdef SPEX_DEBUG
                    // the first entry should be the corresponding pivot
                    ASSERT(P_inv[L->v[j]->i[0]] == j);
                    int sgn;
                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[j]->x[0]));
                    ASSERT (sgn != 0);
#endif
                    for (p = 1; p < L->v[j]->nz; p++)
                    {
                        // row index
                        i = P_inv[L->v[j]->i[p]];
                        if (i > jnext) {continue;}

                        // get the last row-wise off-diagonal entries
                        Lr_offdiag[i] = j; // get the column index
                    }
                }

                Uc_jnext_nz = 0; // # of nnz in U(k+1:jnext-1, jnext)
                for (i = k+1; i <= jnext; i++)
                {
                    ASSERT(Q_inv[U->v[i]->i[0]] == i);
                    for (p = 1 ; p < U->v[i]->nz; p++)
                    {
                        j = Q_inv[U->v[i]->i[p]];
                        if (j > jnext) {continue;}
                        else if (j == jnext)
                        {
                            Uci[Uc_jnext_nz] = i;
                            Ucx[Uc_jnext_nz] = p;
                            Uc_jnext_nz++;
                        }
                        Uc_offdiag[j] = i;
                    }
                }

                LU_scaned = jnext;
            }

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
                    //             x 0 0 0 x 0
                    //             0 x . . . 0
                    //             0 . x . . 0
                    //             0 . . x . 0
                    //             . 0 0 0 x x
                    //                       ^
                    //                       |
                    //                   col n or vk
                    //
                    if (vk_2ndlastnz == -1/*then vk[P[n-1]] must != 0*/  &&
                        Lr_offdiag[n-1] <= k && inext >= n-1)
                    {
                        // perform diagnal swapping with columns n
                        ks = n;
                        SPEX_CHECK(spex_update_finalize_and_insert_vk(vk_dense,
                            h_for_vk, U, L, (const mpz_t*)sd, Q,
                            P_inv, k, n-1, one));
                        SPEX_CHECK(spex_update_dppu1(L, U, sd, Lk_dense_col,
                            Uk_dense_row, &inext, h, Q, Q_inv, P,
                            P_inv, k, ks));
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
                    SPEX_CHECK(spex_update_finalize_and_insert_vk(vk_dense,
                        h_for_vk, U, L, (const mpz_t*)sd, Q,
                        P_inv, k, k, one));
                    SPEX_CHECK(spex_update_cppu(L, U, sd, Lk_dense_col,
                        Uk_dense_row, &inext, &jnext, h, Q, Q_inv, P, P_inv,
                        NULL, NULL, 0, k, ks));
                    break;
                }
            }

            // Since Uc_offdiag is obtained before the updating process, and
            // therefore does not include any potential fillins added to
            // row k of U during this process. For this reason, even if the
            // index of 1st offdiag in column jnext of U is less than k, we
            // will use cppu to swap k and jnext, since we have checked
            // Uk_dense_row->x[Q[jnext]] != 0
            if (inext == k+1 || jnext == k+1 || Uc_offdiag[jnext] <= k)
            {
                if (jnext == k+1 || Uc_offdiag[jnext] <= k)
                {
                    Uc_jnext_nz = 0;
                }
                // use cppu if we see one of the following patterns
                // x . . . .            x 0 0 0 x
                // x x . . .            . x . . 0
                // . . x . .     or     . . x . 0
                // . . . x .            . . . x 0
                // . . . . x            . . . . x
                // These implicitly include the case of jnext == k+1.
                // jnext < n holds since jnext == n has been handled,
                ks = jnext;
                SPEX_CHECK(spex_update_cppu(L, U, sd, Lk_dense_col,
                    Uk_dense_row, &inext, &jnext, h, Q, Q_inv, P, P_inv, Uci,
                    Ucx, Uc_jnext_nz, k, ks));
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
                        maximum = SPEX_MAX(Lr_offdiag[j], Uc_offdiag[j]);
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
                SPEX_CHECK(spex_update_dppu1(L, U, sd, Lk_dense_col,
                    Uk_dense_row, &inext, h, Q, Q_inv, P, P_inv, k, ks));
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
        SPEX_CHECK(spex_update_triangular_solve(vk_dense, h_for_vk, 
            &last_update, NULL /*&vk_2ndlastnz*/, k, L, U,
            (const mpz_t*)sd, P, P_inv));
        // check again in case k is initially n-1
        SPEX_CHECK(SPEX_mpz_sgn(&sgn_vkn, vk_dense->x[P[n-1]]));
        if (sgn_vkn == 0)
        {
            SPEX_FREE_ALL;
            return SPEX_SINGULAR;
        }
        SPEX_CHECK(spex_update_finalize_and_insert_vk(vk_dense, h_for_vk, U, L,
            (const mpz_t*)sd, Q, P_inv, k, k, one));
        // sd[n-1]       = L(P(n-1),n-1)
        SPEX_CHECK(SPEX_mpz_set(sd[n-1],         L->v[n-1]->x[0]));
        // U(n-1,Q(n-1)) = L(P(n-1),n-1)
        SPEX_CHECK(SPEX_mpz_set(U->v[n-1]->x[0], L->v[n-1]->x[0]));
        U->v[n-1]->i[0] = Q[n-1];
        U->v[n-1]->nz = 1;
    }
    else
    {
        // insert U(k,Q(k)) as the first entry in U->v[k]
        p = U->v[k]->nz;
        if (p != 0)
        {
            if (U->v[k]->nzmax <= p)
            {
                // realloc one more entry for U(k,Q[k])
                SPEX_CHECK(SPEX_vector_realloc(U->v[k], U->v[k]->nzmax+1));
            }
            // append U(k,Q[k]) to the end of U->v[k] and swap with the
            // first entry
            SPEX_CHECK(SPEX_mpz_set(U->v[k]->x[p], sd[k]));
            SPEX_CHECK(SPEX_mpz_swap(U->v[k]->x[p], U->v[k]->x[0]));
            U->v[k]->i[p] = U->v[k]->i[0];
            U->v[k]->i[0] = Q[k];
            U->v[k]->nz = p+1;
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_set(U->v[k]->x[0], sd[k]));
            U->v[k]->i[0] = Q[k];
            U->v[k]->nz = 1;
        }
    }
    // S(:,k)=[1;1]
    SPEX_CHECK(SPEX_mpq_set_ui(SL(k), 1, 1));
    SPEX_CHECK(SPEX_mpq_set_ui(SU(k), 1, 1));

#ifdef SPEX_DEBUG
    // check if A=LD^(-1)U
    bool result;
    SPEX_CHECK(spex_update_verify(&result, L, U, A, h, (const mpz_t*) sd,
        P, Q_inv, option));
    printf("the factorization is %s\n", result?"correct":"incorrect");
    if (!result)
    {
        SPEX_FREE_ALL;
        return SPEX_PANIC;
    }
#endif

    SPEX_FREE_ALL;
    return SPEX_OK;
    //return result?SPEX_OK:SPEX_INCORRECT_INPUT;
}
