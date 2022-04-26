//------------------------------------------------------------------------------
// SPEX_Update/SPEX_Update_LU_ColRep.c: perform LU update for column replacement
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis, Erick
// Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

// Purpose: This function performs LU update for column replacement. The
// matrices in the input factorization can be any type and/or kind and does not
// have to be in updatable format. The function will always first check if the
// factorization is updatable and perform necessary conversion if needed. L and
// U in the output factorization will be updatable. The factorization will be
// modified during the update process.  Therefore, if this function fails for
// any reason, the returned F should be considered as undefined.
//
// The matrix A is not modified during the update. If the updated A is
// needed, user can use the follow code if A is in SPEX_dynamic_CSC form.
//
//       SPEX_vector *Vtmp = A->v[k];
//       A->v[k] = vk->v[0];
//       vk->v[0] = Vtmp;

#define SPEX_FREE_ALL                \
    SPEX_FREE(h);                    \
    SPEX_FREE(h_for_vk);             \
    SPEX_FREE(map);                  \
    SPEX_FREE(Lr_offdiag);           \
    SPEX_FREE(Uc_offdiag);           \
    SPEX_FREE(Uci);                  \
    SPEX_FREE(Ucx);                  \
    spex_scattered_vector_free(&Lk_dense_col, option);\
    spex_scattered_vector_free(&Uk_dense_row, option);\
    spex_scattered_vector_free(&vk_dense, option);

#include "spex_update_internal.h"

#define SL(k) (L->v[(k)]->scale)
#define SU(k) (UT->v[(k)]->scale)

SPEX_info SPEX_Update_LU_ColRep
(
    SPEX_factorization* F,  // The SPEX factorization of A, including L, U,
                            // rhos, P, Pinv, Q and Qinv. The factorization
                            // will be modified during the update process.
                            // Therefore, if this function fails for any
                            // reason, the returned F should be considered as
                            // undefined.
    SPEX_matrix *vk,        // Pointer to a n-by-1 dynamic_CSC matrix
                            // which contains the column to be inserted.
                            // The rows of vk are in the same order as A.
    int64_t k,              // The column index that vk will be inserted, 0<=k<n
    const SPEX_options *option// Command parameters
)
{
    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    SPEX_info info;
    if (!spex_initialized()) {return SPEX_PANIC;}

    SPEX_REQUIRE(vk, SPEX_DYNAMIC_CSC, SPEX_MPZ);

    if (!F || F->kind != SPEX_LU_FACTORIZATION || k < 0 || k >= F->L->m ||
        vk->m != F->L->m || vk->n != 1 || vk->v[0]->nz < 0 ||
        (vk->v[0]->nz > 0 && (!(vk->v[0]->x) || !(vk->v[0]->i))))
    {
        return SPEX_INCORRECT_INPUT;
    }
    if (vk->v[0]->nz == 0)
    {
        return SPEX_SINGULAR;
    }

    // make sure F is updatable
    if (!(F->updatable))
    {
        info = SPEX_factorization_convert(F, option);
        if (info != SPEX_OK) return info;
    }

    //--------------------------------------------------------------------------
    // initialize workspace
    //--------------------------------------------------------------------------
    int sgn_vkk, sgn_vkn, r;
    SPEX_matrix *L = F->L, *UT = F->U, *rhos = F->rhos;
    int64_t ks, p, i, j, inext, jnext, n = L->n;
    int64_t *h = NULL, *h_for_vk = NULL, *Lr_offdiag = NULL, *Uc_offdiag = NULL,
        *Uci = NULL, *Ucx = NULL, *map = NULL;
    int64_t *P = F->P_perm, *P_inv = F->Pinv_perm, *Q = F->Q_perm,
        *Q_inv = F->Qinv_perm;
    spex_scattered_vector *Lk_dense_col = NULL, *Uk_dense_row = NULL,
        *vk_dense = NULL;
    mpz_t *sd = rhos->x.mpz;

    //--------------------------------------------------------------------------
    // allocate space for workspace
    //--------------------------------------------------------------------------
    h          = (int64_t*) SPEX_calloc(n, sizeof(int64_t));
    h_for_vk   = (int64_t*) SPEX_calloc(n, sizeof(int64_t));
    Uc_offdiag = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    Lr_offdiag = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    Uci        = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    Ucx        = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    map        = (int64_t*) SPEX_malloc(n* sizeof(int64_t));
    if (!h || !h_for_vk || !Uci || !Ucx || !map || !Uc_offdiag || !Lr_offdiag)
    {
        SPEX_FREE_ALL;
        return SPEX_OUT_OF_MEMORY;
    }

    //--------------------------------------------------------------------------
    // initialize environment for the inserted column
    //--------------------------------------------------------------------------
    SPEX_CHECK(spex_update_get_scattered_v(&vk_dense, NULL, vk->v[0], n, n,
        NULL, true, option));
    //index of column of L used to perform IPGE for vk
    int64_t last_update = -1;
    // the index of last nnz in vk(P[last_update+1,n-2])
    int64_t vk_2ndlastnz = -1;
    // P_inv[vk->i[0...vk_top-1]] <= last_update
    int64_t vk_top = 0;
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

    //--------------------------------------------------------------------------
    // initialize for the while loop
    //--------------------------------------------------------------------------
    LU_scaned = Q_inv[k];
    // remove column k of U
    for (i = 0; i < LU_scaned; i++)
    {
        for (p = 1 ; p < UT->v[i]->nz; p++)
        {
            j = UT->v[i]->i[p];
            if (j == k)
            {
                // move the last entry to current position
                UT->v[i]->nz--;
                SPEX_CHECK(SPEX_mpz_swap(UT->v[i]->x[p],
                                         UT->v[i]->x[UT->v[i]->nz]));
                UT->v[i]->i[p] = UT->v[i]->i[UT->v[i]->nz];
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
        // build Lk_dense_col and Uk_dense_row, remove explicit 0, search for
        // inext and jnext
        SPEX_CHECK(spex_update_get_scattered_v(&Lk_dense_col, &inext,
            L->v[k], n, k, P_inv, false, option));
        SPEX_CHECK(spex_update_get_scattered_v(&Uk_dense_row, &jnext,
            UT->v[k], n, k, Q_inv, false, option));
    }

    //--------------------------------------------------------------------------
    // push column k to position n-1
    //--------------------------------------------------------------------------
    while (k < n-1)
    {
#ifdef SPEX_DEBUG
        // check for L and Lk_dense_col
        for (int64_t ii = 0; ii < n; ii++)
        {
            // check if the first entry of all L->v is the pivot
            if (ii != k && L->v[ii]->i[0] != P[ii])
            {
                printf("Incorrect col %ld of L\n",ii);
                SPEX_CHECK(SPEX_PANIC);
            }

            // check if Lk_dense_col has nonzero entry but not in nnz pattern
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

        // check for UT and Uk_dense_row
        for (int64_t ii =0;ii<n;ii++)
        {
            // check if the first entry of all UT->v is the pivot
            if (ii != k && UT->v[ii]->i[0] != Q[ii])
            {
                printf("Incorrect row %ld of U\n",ii);
                SPEX_CHECK(SPEX_PANIC);
            }

            // check if Uk_dense_row has nonzero entry but not in nnz pattern
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

        // check for sd and h
        for (int64_t ii = 0; ii < n; ii++)
        {
            // check if any pivot is zero
            SPEX_CHECK(SPEX_mpz_sgn(&r, sd[ii]));
            if (r == 0)
            {
                printf("sd[%ld]=0 at file %s line %d\n",ii,__FILE__,__LINE__);
                SPEX_CHECK(SPEX_PANIC);
            }

            // check if any entry in history vector is out of range
            if (h[ii] < -1)
            {
                printf("h[%ld]=%ld at file %s line %d\n", ii, h[ii],
                    __FILE__, __LINE__);
                SPEX_CHECK(SPEX_PANIC);
            }
        }

        if (jnext < n)
        {
            ASSERT (jnext > k);
            SPEX_CHECK(SPEX_mpz_sgn(&r, Uk_dense_row->x[Q[jnext]]));
            ASSERT (r != 0);
        }
        else // jnext == n, i.e., only 1 nnz in row k of U
        {
            // explicit zero is always removed from nnz pattern of Uk_dense_row
            // but not for Lk_dense_col
            ASSERT(Uk_dense_row->nz == 1);
        }
        if (inext < n)
        {
            ASSERT (inext > k);
            SPEX_CHECK(SPEX_mpz_sgn(&r, Lk_dense_col->x[P[inext]]));
            ASSERT (r != 0);
        }
#endif

        //----------------------------------------------------------------------
        // get the (k-1)-th IPGE update of inserted column, so that vk
        // will be ready to be inserted as column k if needed
        //----------------------------------------------------------------------
        SPEX_CHECK(spex_update_triangular_solve(vk_dense, &vk_top, h_for_vk,
            &last_update, &vk_2ndlastnz, k, L, UT,
            (const SPEX_matrix*)rhos, P, P_inv));
        SPEX_CHECK(SPEX_mpz_sgn(&sgn_vkk, vk_dense->x[P[k]]));
        SPEX_CHECK(SPEX_mpz_sgn(&sgn_vkn, vk_dense->x[P[n-1]]));

        //----------------------------------------------------------------------
        // check for singularity
        //----------------------------------------------------------------------
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
        // case 1
        if (jnext == n && sgn_vkk == 0 )
        {
            SPEX_FREE_ALL;
            return SPEX_SINGULAR;
        }
        // case 2
        if (vk_2ndlastnz == -1 && sgn_vkn == 0)
        {
            SPEX_FREE_ALL;
            return SPEX_SINGULAR;
        }
        // case 3, check only when whole U is scaned 
        if (vk_2ndlastnz == -1 && LU_scaned == n-1 && Uc_offdiag[n-1] <= k)
        {
            SPEX_CHECK(SPEX_mpz_sgn(&r, Uk_dense_row->x[Q[n-1]]));
            if (r == 0)
            {
                SPEX_FREE_ALL;
                return SPEX_SINGULAR;
            }
        }

        //----------------------------------------------------------------------
        // see if we can finish the update here with column swapping
        //----------------------------------------------------------------------
        // if the next nnz in current row is in vk, then use vk, which will
        // help to avoid performing extra ipge iterations for vk
        if (jnext == n)
        {
            // if jnext == n, swapping columns k and vk will be more efficient,
            // since there is no need to backtrack vk (i.e., col n) (we just
            // perform k-th IPGE iteration for column n), and column n-1 can be
            // updated by scaling after CPPU.
            ks = n;
            SPEX_CHECK(spex_update_finalize_and_insert_vk(vk_dense, h_for_vk,
                UT, L, (const SPEX_matrix*)rhos, Q, P_inv, k, k, option));

            SPEX_CHECK(spex_update_cppu(L, UT, rhos, Lk_dense_col,
                Uk_dense_row, &inext, &jnext, h, Q, Q_inv,
                P, P_inv, NULL, NULL, 0, k, ks, option));
            break;
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
                    ASSERT(Q_inv[UT->v[i]->i[0]] == i);
                    for (p = 1 ; p < UT->v[i]->nz; p++)
                    {
                        j = Q_inv[UT->v[i]->i[p]];
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
                    h_for_vk, UT, L, (const SPEX_matrix*)rhos, Q, P_inv, k,
                    n-1, option));
                ks = n;
            }
            if (jnext > ks || (ks == n && jnext >= n-1 && sgn_vkk == 0))
            {
                SPEX_CHECK(spex_update_dppu1(L, UT, rhos, Lk_dense_col,
                    Uk_dense_row, &inext, h, Q, Q_inv, P, P_inv, k, ks,
                    option));
            }
            else
            {
                SPEX_CHECK(spex_update_dppu2(L, UT, rhos, Lk_dense_col,
                    Uk_dense_row, &jnext, h, Q, Q_inv, P, P_inv, k, ks,
                    option));
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
                    ASSERT(Q_inv[UT->v[i]->i[0]] == i);
                    for (p = 1 ; p < UT->v[i]->nz; p++)
                    {
                        j = Q_inv[UT->v[i]->i[p]];
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
                            h_for_vk, UT, L, (const SPEX_matrix*)rhos, Q,
                            P_inv, k, n-1, option));
                        SPEX_CHECK(spex_update_dppu1(L, UT, rhos, Lk_dense_col,
                            Uk_dense_row, &inext, h, Q, Q_inv, P,
                            P_inv, k, ks, option));
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
                        h_for_vk, UT, L, (const SPEX_matrix*)rhos, Q,
                        P_inv, k, k, option));
                    SPEX_CHECK(spex_update_cppu(L, UT, rhos, Lk_dense_col,
                        Uk_dense_row, &inext, &jnext, h, Q, Q_inv, P, P_inv,
                        NULL, NULL, 0, k, ks, option));
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
                SPEX_CHECK(spex_update_cppu(L, UT, rhos, Lk_dense_col,
                    Uk_dense_row, &inext, &jnext, h, Q, Q_inv, P, P_inv, Uci,
                    Ucx, Uc_jnext_nz, k, ks, option));
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
                SPEX_CHECK(spex_update_dppu1(L, UT, rhos, Lk_dense_col,
                    Uk_dense_row, &inext, h, Q, Q_inv, P, P_inv, k, ks,
                    option));
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

    //--------------------------------------------------------------------------
    // perform any remaining work to finish the update
    //--------------------------------------------------------------------------
    // k will be the column index where vk is inserted
    if (k == n-1)
    {
        SPEX_CHECK(spex_update_triangular_solve(vk_dense, &vk_top, h_for_vk, 
            &last_update, NULL /*&vk_2ndlastnz*/, k, L, UT,
            (const SPEX_matrix*)rhos, P, P_inv));
        // check again in case k is initially n-1
        SPEX_CHECK(SPEX_mpz_sgn(&sgn_vkn, vk_dense->x[P[n-1]]));
        if (sgn_vkn == 0)
        {
            SPEX_FREE_ALL;
            return SPEX_SINGULAR;
        }
        SPEX_CHECK(spex_update_finalize_and_insert_vk(vk_dense, h_for_vk, UT, L,
            (const SPEX_matrix*)rhos, Q, P_inv, k, k, option));
        // sd[n-1]       = L(P(n-1),n-1)
        SPEX_CHECK(SPEX_mpz_set(sd[n-1],         L->v[n-1]->x[0]));
        // U(n-1,Q(n-1)) = L(P(n-1),n-1)
        SPEX_CHECK(SPEX_mpz_set(UT->v[n-1]->x[0], L->v[n-1]->x[0]));
        UT->v[n-1]->i[0] = Q[n-1];
        UT->v[n-1]->nz = 1;
    }
    else
    {
        // insert U(k,Q(k)) as the first entry in U->v[k]
        p = UT->v[k]->nz;
        if (p != 0)
        {
            if (UT->v[k]->nzmax <= p)
            {
                // realloc one more entry for U(k,Q[k])
                SPEX_CHECK(SPEX_vector_realloc(UT->v[k], UT->v[k]->nzmax+1,
                    option));
            }
            // append U(k,Q[k]) to the end of U->v[k] and swap with the
            // first entry
            SPEX_CHECK(SPEX_mpz_set(UT->v[k]->x[p], sd[k]));
            SPEX_CHECK(SPEX_mpz_swap(UT->v[k]->x[p], UT->v[k]->x[0]));
            UT->v[k]->i[p] = UT->v[k]->i[0];
            UT->v[k]->i[0] = Q[k];
            UT->v[k]->nz = p+1;
        }
        else
        {
            SPEX_CHECK(SPEX_mpz_set(UT->v[k]->x[0], sd[k]));
            UT->v[k]->i[0] = Q[k];
            UT->v[k]->nz = 1;
        }
    }
    // S(:,k)=[1;1]
    SPEX_CHECK(SPEX_mpq_set_ui(SL(k), 1, 1));
    SPEX_CHECK(SPEX_mpq_set_ui(SU(k), 1, 1));

    SPEX_FREE_ALL;
    return SPEX_OK;
}
