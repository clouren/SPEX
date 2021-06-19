//------------------------------------------------------------------------------
// SPEX_Update/Test/ptest.c: performance test for SPEX_LU_Update library
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_Update/License for the license.

//------------------------------------------------------------------------------

/*
 * performance test
 */

#define FREE_WORKSPACE                           \
{                                                \
    fclose(out_file);                            \
    SPEX_FREE(option);                           \
    SPEX_finalize() ;                            \
}

#include "test.h"
#include <assert.h>
int64_t init_basis[27]={27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 62, 61, 66, 46, 50, 44, 45, 48, 49, 59, 60, 64, 65, 76, 75};


int main( int argc, char* argv[])
{
    //char *prob_name = "lp_80bau3b";
    char *prob_name = "lp_25fv47";//5.5018458883E+03
    //char *prob_name = "lp_afiro"; // optimal: -4.6475314E+02
    //char *prob_name = "aa5";
    if (argc >= 2)
    {
        prob_name = argv[1];
    }
    SPEX_info info;
    //------------------------------------------------------------------
    // Initialize for SPEX libray
    //------------------------------------------------------------------

    SPEX_initialize () ;

    //------------------------------------------------------------------
    // Allocate memory
    //------------------------------------------------------------------

    int64_t n, i, j, p, nz;
    int sgn;
    double t1 = 0, t2 = 0, t3 = 0, t_solve = 0;
    SPEX_options* option = NULL;
    SPEX_matrix *Prob_A = NULL, *Prob_c = NULL, *tempA = NULL;
    SPEX_matrix *L1 = NULL, *U1 = NULL, *rhos = NULL, *A1 = NULL;
    SPEX_matrix *rhos2 = NULL, *rhos3 = NULL;
    SPEX_matrix *L2 = NULL, *U2 = NULL, *A2 = NULL;
    SPEX_matrix *b = NULL, *c = NULL,  *x  = NULL, *y = NULL;
    SPEX_matrix *L3 = NULL, *U3 = NULL, *A3 = NULL;
    SPEX_vector *tmpv, *vk = NULL, *vk3 = NULL;
    mpz_t z, minz, tmpz;
    mpq_t minq, tmpq, one;
    SPEX_LU_analysis* analysis = NULL;
    int64_t *P1_inv = NULL, *P = NULL, *P_inv = NULL, *Q = NULL, *Q_inv = NULL;
    int64_t *P3 = NULL, *P3_inv = NULL, *Q3 = NULL, *Q3_inv = NULL;
    int64_t *basis = NULL, *h = NULL, *used_as_basis = NULL;
    clock_t start_llu, start_luu, end_llu, end_luu, start1, end1, start2, end2,
            start_luu3 = 0, end_luu3 = 0;
    FILE *out_file = NULL;
    glp_prob *LP;
    glp_smcp parm;
    LP = glp_create_prob();
    glp_init_smcp(&parm);
    parm.it_lim = 1;
    OK(SPEX_create_default_options(&option));
    OK(SPEX_mpz_init(z));
    OK(SPEX_mpz_init(minz));
    OK(SPEX_mpz_init(tmpz));
    OK(SPEX_mpq_init(minq));
    OK(SPEX_mpq_init(tmpq));
    OK(SPEX_mpq_init(one));
    OK(SPEX_mpq_set_ui(one, 1, 1));

    //--------------------------------------------------------------------------
    // open output file
    //--------------------------------------------------------------------------
    out_file = fopen("result.out", "w");
    if (out_file == NULL)
    {
        printf("Could not open file");
        return 0;
    }

    //--------------------------------------------------------------------------
    // read in matrix
    //--------------------------------------------------------------------------
    char file_name[1000] = "TestMats/";
    strcat(file_name, prob_name);
    // read matrix A, arrays b, lb, ub and c to LP
    SPEX_construct_LP(LP, &Prob_A, &b, &Prob_c, file_name, option);
    n = Prob_A->m;

    //--------------------------------------------------------------------------
    // generate initial basis matrix
    //--------------------------------------------------------------------------
    glp_adv_basis(LP, 0);
    glp_factorize(LP); // in order to use glp_get_bhead
    basis = (int64_t*) SPEX_malloc(n * sizeof(int64_t));
    used_as_basis = (int64_t*) SPEX_malloc(Prob_A->n * sizeof(int64_t));
    if (!basis || !used_as_basis)
    {
        FREE_WORKSPACE;
        return 0;
    }
    for (i = 0; i < Prob_A->n; i++)
    {
        used_as_basis[i] = -1;
    }

    printf("getting initial basic variables....\n");
    parm.msg_lev = GLP_MSG_OFF;
    while (glp_get_status(LP) != GLP_FEAS)
    {
        glp_simplex(LP, &parm);
    }
    for (i = 0; i < n; i++)
    {
#if 0
        basis[i] = init_basis[i];
#else
        basis[i] = glp_get_bhead(LP, i+1)-1;
        if (basis[i] < n && glp_get_mat_row(LP, basis[i]+1, NULL, NULL) != 0)
        {
            //printf("basis[%ld]=%ld\n",i,basis[i]);
            GOTCHA;
            return 0;;
        }
#endif
    }
#if 0
    //while (glp_get_status(LP) != GLP_OPT)
    //{
    //    glp_simplex(LP, &parm);
    double opt_val =0.0;
    for (j = 0; j<Prob_A->n;j++)
    {
        if (glp_get_col_stat(LP,j+1)==GLP_BS)
        opt_val += glp_get_col_prim(LP,j+1)*glp_get_obj_coef(LP,j+1);
        if (glp_get_col_prim(LP,j+1)<0)
            printf("panic %ld %f\n",j, glp_get_col_prim(LP,j+1));
    }
    //printf("%f\n",opt_val);
    //}
    //return 0;
#endif

    OK(SPEX_mpq_get_den(tmpz, Prob_A->scale));
    OK(SPEX_mpz_cmp_ui(&sgn, tmpz, 1));
    if (sgn != 0)
    {
        OK(SPEX_gmp_printf("scale is %Qd, whose den is not 1\n",
            Prob_A->scale));
        FREE_WORKSPACE;
        return 0;
    }

    // allocate A2 with n sparse vectors with initially 0 nnz
    OK(SPEX_matrix_allocate(&A2, SPEX_DYNAMIC_CSC, SPEX_MPZ, n, n, 0,
        false, true, option));
    for (i = 0; i < n; i++)
    {
        if (basis[i] < n)
        {
            OK(SPEX_vector_realloc(A2->v[i], 1, option));
            mpq_get_num(A2->v[i]->x[0], Prob_A->scale);
            A2->v[i]->i[0] = basis[i];
            A2->v[i]->nz = 1;
        }
        else
        {
            j = basis[i]-n;
            used_as_basis[j] = i;
            OK(SPEX_vector_realloc(A2->v[i], Prob_A->p[j+1]-Prob_A->p[j],
                option));
            nz = 0;
            for (p = Prob_A->p[j]; p < Prob_A->p[j+1]; p++)
            {
                A2->v[i]->i[nz] = Prob_A->i[p];
                OK(SPEX_mpz_set(A2->v[i]->x[nz], SPEX_1D(Prob_A, p, mpz)));
                nz++;
            }
            A2->v[i]->nz = nz;
        }
    }
    OK(SPEX_mpq_set(A2->scale, Prob_A->scale));
    // get a CSC x MPZ A1 from dynamic_CSC x MPZ A2
    OK(SPEX_matrix_copy(&A1, SPEX_CSC, SPEX_MPZ, A2, option));

    //--------------------------------------------------------------------------
    // perform LU factorization for the initial matrix A1
    //--------------------------------------------------------------------------
    start_llu = clock();

    // We now perform symbolic analysis by getting the column preordering of A
    OK(SPEX_LU_analyze(&analysis, A1, option));

    // Now we perform the SPEX Left LU factorization to obtain matrices L and U
    // and a row permutation P such that PAQ = LDU.
    OK(SPEX_Left_LU_factorize(&L1, &U1, &rhos, &P1_inv, A1, analysis, option));

    end_llu = clock();

    //--------------------------------------------------------------------------
    // generate initial inputs for LU update
    //--------------------------------------------------------------------------
    // generate permutation vectors P, Q, P_inv, Q_inv and vectors rhos
    Q      = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    Q3     = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    P      = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    P3     = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    Q_inv  = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    Q3_inv = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    P_inv  = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    P3_inv = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    h  = (int64_t*) SPEX_calloc(n, sizeof(int64_t));
    if (!P  || !Q  || !P_inv  || !Q_inv  ||
        !P3 || !Q3 || !P3_inv || !Q3_inv || !h)
    {
        FREE_WORKSPACE;
        return 0;
    }
    OK(SPEX_matrix_allocate(&rhos2, SPEX_DENSE, SPEX_MPZ, n, 1, n, false, true,
        option));
    OK(SPEX_matrix_allocate(&rhos3, SPEX_DENSE, SPEX_MPZ, n, 1, n, false, true,
        option));
    for (i = 0; i < n; i++)
    {
        Q[i] = analysis->q[i];
        P_inv[i] = P1_inv[i];
        P[P_inv[i]] = i;
        Q_inv[Q[i]] = i;
        OK(SPEX_mpz_set(SPEX_1D(rhos2, i, mpz), SPEX_1D(rhos, i, mpz)));
    }

    // convert the factorization to SPEX_mat to be used in the update process
    //OK(SPEX_CSC_to_mat(&L2, P, true,  L1, option));
    OK(SPEX_matrix_copy(&L2, SPEX_DYNAMIC_CSC, SPEX_MPZ, L1, option));
    OK(SPEX_Update_permute_row(L2, P, option));
    OK(SPEX_Update_matrix_canonicalize(L2, P, option));
    //OK(SPEX_CSC_to_mat(&U2, Q, false, U1, option));// U2 stored in row-wise
    OK(SPEX_transpose(&tempA, U1, option));
    OK(SPEX_matrix_copy(&U2, SPEX_DYNAMIC_CSC, SPEX_MPZ, tempA, option));
    OK(SPEX_Update_permute_row(U2, Q, option));
    OK(SPEX_Update_matrix_canonicalize(U2, Q, option));

    // allocate space for vk
    OK(SPEX_vector_allocate(&vk, 0, option));
    OK(SPEX_vector_allocate(&vk3, 0, option));

    // allocate space for c and y
    OK(SPEX_matrix_allocate(&c, SPEX_DENSE, SPEX_MPZ, n, 1, n, false, true,
        option));
    OK(SPEX_matrix_allocate(&y, SPEX_DENSE, SPEX_MPZ, n, 1, n, false, true,
        option));

    printf("scales of A and b are %s\n", mpq_equal(Prob_A->scale, b->scale)==0?
        "different":"same");
    int64_t k = 0, new_col = 0, iter = 0;
    while (1)
    {
        int64_t L_nnz = 0, U_nnz = 0, L3_nnz = 0, U3_nnz = 0;
        size_t U1_sum_size = 0, U2_sum_size = 0, U3_sum_size = 0,
               L1_sum_size = 0, L2_sum_size = 0, L3_sum_size = 0, s;
        t_solve = 0;
        start1 = clock();
        start2 = clock();
        // solve for x_basic
        OK(SPEX_Update_Solve(&x, b, L2, U2, A2->scale, h,
            (const SPEX_matrix*)rhos2, P, Q_inv, option));
        end2 = clock();
        t_solve += (double) (end2 - start2) / CLOCKS_PER_SEC;

        // reset objective value z = 0
        OK(SPEX_mpz_set_ui(z, 0));
        for (i = 0; i < n; i++)
        {
            j = basis[i];
            if (j < n)
            {
                //j = j-1;
                // build vector c
                OK(SPEX_mpz_set_ui(c->x.mpz[i], 0));
            }
            else
            {
                j = j-n;
                // build vector c
                OK(SPEX_mpz_set(c->x.mpz[i], Prob_c->x.mpz[j]));

                // compute objective value z
                OK(SPEX_mpz_addmul(z, c->x.mpz[i], x->x.mpz[i]));

                OK(SPEX_mpq_set_z(tmpq, x->x.mpz[i]));
                OK(SPEX_mpq_div(tmpq, tmpq, x->scale));
#if 1
                if (mpz_sgn(x->x.mpz[i]) * mpq_sgn(x->scale) < 0)
                {
                    printf("xz[%ld]<0 %d %d %f \n",j,mpz_sgn(x->x.mpz[i]),mpq_sgn(x->scale),mpq_get_d(tmpq));
                    gmp_fprintf(out_file,"%Qd\n",tmpq);
                    OK(SPEX_PANIC);
                }
#else
                SPEX_gmp_printf("%ld xz[%ld]= %Qd \n",i,j,tmpq);
#endif
            }
        }

        // write current basis to output file
        for (i = 0 ; i < n; i++)
        {
            fprintf(out_file, "%ld, ",basis[i]);
        }
        fprintf(out_file,"\n");

        // set c->scale = Prob_c->scale
        OK(SPEX_mpq_set(c->scale, Prob_c->scale));

        // compute the real objective value with scales applied
        OK(SPEX_mpq_set_z(tmpq, z));
        OK(SPEX_mpq_div(tmpq, tmpq, x->scale));
        OK(SPEX_mpq_div(tmpq, tmpq, Prob_c->scale));
        //OK(SPEX_gmp_printf("obj value = %Qd\n",tmpq));
        printf("obj value = %f\n", mpq_get_d(tmpq));

        //----------------------------------------------------------------------
        // find the entering variable
        //----------------------------------------------------------------------
        // solve A'*c_new = c for updated coefficient for objective function
        start2 = clock();
        OK(SPEX_Update_Solve(&c, c, U2, L2, A2->scale, h,
            (const SPEX_matrix*) rhos2, Q, P_inv, option));
        end2 = clock();
        t_solve += (double) (end2 - start2) / CLOCKS_PER_SEC;

        new_col = -1;
        OK(SPEX_mpz_set_ui(minz, 0));// find the most negative coefficient
        double mind = 0.0, tmpd = 0.0;
        for (j = 0; j < Prob_A->n; j++)
        {
            if (used_as_basis[j] >= 0) {continue;} // skip basis variables

            // compute tmpz = -c[j]+c_new'*A(:,j)
            OK(SPEX_mpz_set_ui(tmpz, 0));
            // iterate across nnz in column j of A
            for (p = Prob_A->p[j]; p < Prob_A->p[j+1]; p++)
            {
                i = Prob_A->i[p];
                OK(SPEX_mpz_addmul(tmpz, Prob_A->x.mpz[p], c->x.mpz[i]));
            }
            OK(SPEX_mpq_set_z(tmpq, tmpz));
            OK(SPEX_mpq_div(tmpq, tmpq, Prob_A->scale));
            OK(SPEX_mpq_div(tmpq, tmpq,      c->scale));
            tmpd = mpq_get_d(tmpq);
            OK(SPEX_mpq_set_z(tmpq, Prob_c->x.mpz[j]));
            OK(SPEX_mpq_div(tmpq, tmpq, Prob_c->scale));
            tmpd = tmpd - mpq_get_d(tmpq);
            //OK(SPEX_mpz_cmp(&sgn, tmpz, minz));
            //if (sgn < 0)
            if (tmpd > mind)//TODO
            {
                //OK(SPEX_gmp_printf("c:%Zd <%Zd\n",tmpz,minz));
                //OK(SPEX_mpz_set(minz, tmpz));
                //printf("c: %f < %f new_col = %ld --> %ld\n",tmpd, mind,new_col,j);
                mind = tmpd;
                new_col = j;
            }
        }

        if (new_col == -1)
        {
            printf("optimal solution found!\n");
            break;
        }

        //----------------------------------------------------------------------
        // construct vk and find the existing variable
        //----------------------------------------------------------------------
        // reset y[i] = 0
        for (i = 0; i < n; i++)
        {
            OK(SPEX_mpz_set_ui(y->x.mpz[i], 0));
        }
        // allocate space for vk if needed
        vk->nz = Prob_A->p[new_col+1]-Prob_A->p[new_col];
        if (vk->nzmax < vk->nz)
        {
            OK(SPEX_vector_realloc(vk, vk->nz, option));
        }
        if (vk3->nzmax < vk->nz)
        {
            OK(SPEX_vector_realloc(vk3, vk->nz, option));
        }
        i = 0;
        j = 0;
        // iterate across nnz in the entering column of A
        for (p = Prob_A->p[new_col]; p < Prob_A->p[new_col+1]; p++)
        {
            for (; j < Prob_A->i[p]; j++)
            {
                OK(SPEX_mpz_set_ui(y->x.mpz[j], 0));
            }
            OK(SPEX_mpz_set(y->x.mpz[j], SPEX_1D(Prob_A, p, mpz)));// dense y

            vk->i[i] = j;
            OK(SPEX_mpz_set(vk->x[i], SPEX_1D(Prob_A, p, mpz)));// sparse vk
            vk3->i[i] = j;
            OK(SPEX_mpz_set(vk3->x[i], SPEX_1D(Prob_A, p, mpz)));// sparse vk
            i++;
            j++;

        }
        vk->nz = i;
        vk3->nz = i;
        // set y->scale = Prob_A->scale
        OK(SPEX_mpq_set(y->scale, Prob_A->scale));
#if 0
        printf("\n double checking the nnz in column %ld:\n",new_col);
        for (i = 0; i<n;i++)
        {
            if (mpz_sgn(y->x.mpz[i])!=0)printf("%ld ",i);
        }
        printf("\n nnz in solution:\n");
#endif

        // solve for Ay = y
        start2 = clock();
        OK(SPEX_Update_Solve(&y, y, L2, U2, A2->scale, h,
            (const SPEX_matrix*) rhos2, P, Q_inv, option));
        end2 = clock();
        t_solve += (double) (end2 - start2) / CLOCKS_PER_SEC;

        // perform ratio test to find existing variable
        k = -1;
        int sd_sgn = mpz_sgn(rhos2->x.mpz[n-1]);
        for (i = 0; i < n; i++)
        {
            int y_sgn = mpz_sgn(y->x.mpz[i]) * sd_sgn;

            if (y_sgn > 0)
            {
                if (mpz_sgn(x->x.mpz[i]) == 0)
                {
                    k = i;
                    printf("found 0!! This update won't improve obj value\n");
                    break;
                }

                OK(SPEX_mpq_set_num(tmpq, x->x.mpz[i]));
                OK(SPEX_mpq_set_den(tmpq, y->x.mpz[i]));
                OK(SPEX_mpq_canonicalize(tmpq));
                OK(SPEX_mpq_mul(tmpq, tmpq, y->scale));
                OK(SPEX_mpq_div(tmpq, tmpq, x->scale));
                        //printf(" x/y=%f \n",mpq_get_d(tmpq));
                if (k == -1)
                {
                    OK(SPEX_mpq_set(minq, tmpq));
                    k = i;
                }
                else
                {
                    OK(SPEX_mpq_cmp(&sgn, tmpq, minq));
                    if (sgn < 0)
                    {
                        //OK(SPEX_gmp_printf("%Qd <%Qd\n",tmpq,minq));
                        //printf("%f < %f\n",mpq_get_d(tmpq),mpq_get_d(minq));
                        OK(SPEX_mpq_set(minq, tmpq));
                        k = i;
                    }
                }
            }
        }
        if (k == -1)
        {
            printf("LP is unbound!\n");
            return 0;
        }
        end1 = clock();
        fprintf(out_file,"prev basis k(%ld): %ld; new basis k: %ld\n",
            k, basis[k]-n, new_col);
        printf("prev basis k(%ld): %ld; new basis k: %ld\n",
            k, basis[k]-n, new_col);
        used_as_basis[new_col] = k;
        used_as_basis[basis[k]-n] = -1;
        basis[k] = new_col+n;
        fprintf(out_file,
            "\n---------------------------------------------------------\n");
        fprintf(out_file,
            "----------%ld: replacing k(%ld) with new_col(%ld)-----------\n",
            iter, k, new_col);
        fprintf(out_file,
            "-----------------------------------------------------------\n");
        printf("\n---------------------------------------------------------\n");
        printf("----------%ld: replacing k(%ld) with new_col(%ld)-----------\n",
            iter, k, new_col);
        printf("-----------------------------------------------------------\n");
        printf("\nSolving 3 linear equations time: \t%lf", t_solve);
        printf("\nSearch k and new_col time: \t\t%lf\n",
           (double) (end1 - start1) / CLOCKS_PER_SEC);

        //----------------------------------------------------------------------
        // generate new matrix with vk inserted
        //----------------------------------------------------------------------
        tmpv = A2->v[k]; A2->v[k] = vk; vk = tmpv;
        OK(SPEX_matrix_free(&A1, option));
        OK(SPEX_matrix_copy(&A1, SPEX_CSC, SPEX_MPZ, A2, option));
        if (iter == 0)
        {
            OK(SPEX_matrix_copy(&A3, SPEX_DYNAMIC_CSC, SPEX_MPZ, A2, option));
        }
        tmpv = A2->v[k]; A2->v[k] = vk; vk = tmpv;

        //----------------------------------------------------------------------
        // perform LU factorization for matrix A1
        //----------------------------------------------------------------------
        SPEX_matrix_free(&L1, NULL);
        SPEX_matrix_free(&U1, NULL);
        SPEX_matrix_free(&rhos, NULL);
        SPEX_FREE(P1_inv);
        SPEX_FREE(analysis);
        start_llu = clock();

        // perform symbolic analysis by getting the column preordering of A
        OK(SPEX_LU_analyze(&analysis, A1, option));

        // Now we perform the SPEX Left LU factorization to obtain matrices L
        // and U and a row permutation P such that PAQ = LDU.
        OK(SPEX_Left_LU_factorize(&L1, &U1, &rhos, &P1_inv, A1, analysis,
            option));
//        if (info == SPEX_OK) {printf("matrix is not singular!\n");}

        end_llu = clock();

        //----------------------------------------------------------------------
        // perform LU update using the factorization from last loop
        //----------------------------------------------------------------------

        if (iter > 0)
        {
            start_luu3 = clock();

            OK(SPEX_Update_LU_ColRep(A3, L3, U3, rhos3, P3, P3_inv, Q3, Q3_inv,
                &vk3, k, option));

            end_luu3 = clock();

            for (i = 0; i < n; i++)
            {
                L3_nnz += L3->v[i]->nz;
                for (p = 0; p < L3->v[i]->nz; p++)
                {
                    OK(SPEX_mpz_sizeinbase(&s, L3->v[i]->x[p], 2));
                    L3_sum_size += s;
                }
                U3_nnz += U3->v[i]->nz;
                for (p = 0; p < U3->v[i]->nz; p++)
                {
                    OK(SPEX_mpz_sizeinbase(&s, U3->v[i]->x[p], 2));
                    U3_sum_size += s;
                }
            }
        }
        for (i = 0; i < n; i++)
        {
            Q3[i] = analysis->q[i];
            P3_inv[i] = P1_inv[i];
            P3[P3_inv[i]] = i;
            Q3_inv[Q3[i]] = i;
            OK(SPEX_mpz_set(SPEX_1D(rhos3, i, mpz), SPEX_1D(rhos, i, mpz)));
        }
        // convert L and U to SPEX_mat to be used in the update process
        SPEX_matrix_free(&L3, option);
        SPEX_matrix_free(&U3, option);
        //OK(SPEX_CSC_to_mat(&L3, P3, true,  L1, option));
        OK(SPEX_matrix_copy(&L3, SPEX_DYNAMIC_CSC, SPEX_MPZ, L1, option));
        OK(SPEX_Update_permute_row(L3, P3, option));
        OK(SPEX_Update_matrix_canonicalize(L3, P3, option));
        //OK(SPEX_CSC_to_mat(&U3, Q3, false, U1, option));// U2 stored in row-wise
        OK(SPEX_transpose(&tempA, U1, option));
        OK(SPEX_matrix_copy(&U3, SPEX_DYNAMIC_CSC, SPEX_MPZ, tempA, option));
        OK(SPEX_matrix_free(&tempA, option));
        OK(SPEX_Update_permute_row(U3, Q3, option));
        OK(SPEX_Update_matrix_canonicalize(U3, Q3, option));

        //----------------------------------------------------------------------
        // perform LU update for matrix A2->A1
        //----------------------------------------------------------------------
        start_luu = clock();

        OK(SPEX_Update_LU_ColRep(A2, L2, U2, rhos2, P, P_inv, Q, Q_inv, &vk, k,
            option));

        end_luu = clock();

        //----------------------------------------------------------------------
        // print results
        //----------------------------------------------------------------------
        // Timing stats
        double t_llu = (double) (end_llu-start_llu)/CLOCKS_PER_SEC;
        double t_luu3 = (double) (end_luu3 - start_luu3) / CLOCKS_PER_SEC;
        double t_luu = (double) (end_luu - start_luu) / CLOCKS_PER_SEC;
        if (iter > 0)
        {
            t1 += t_llu;
            t2 += t_luu;
            t3 += t_luu3;
        }
        for (i = 0; i < n; i++)
        {
            L_nnz += L2->v[i]->nz;
            for (p = 0; p < L2->v[i]->nz; p++)
            {
                OK(SPEX_mpz_sizeinbase(&s, L2->v[i]->x[p], 2));
                L2_sum_size += s;
            }
            U_nnz += U2->v[i]->nz;
            for (p = 0; p < U2->v[i]->nz; p++)
            {
                OK(SPEX_mpz_sizeinbase(&s, U2->v[i]->x[p], 2));
                U2_sum_size += s;
            }
        }
        for (p = 0; p < L1->p[n]; p++)
        {
            OK(SPEX_mpz_sizeinbase(&s, L1->x.mpz[p], 2));
            L1_sum_size += s;
        }
        for (p = 0; p < U1->p[n]; p++)
        {
            OK(SPEX_mpz_sizeinbase(&s, U1->x.mpz[p], 2));
            U1_sum_size += s;
        }
        printf("\n\ttime     \tL_nnz \tsum(Lb) \tU_nnz \t\tsum(Ub)\n");
        printf("\nLLU: \t%lf \t%ld   \t%ld    \t%ld   \t%ld", t_llu,
            L1->p[n], L1_sum_size, U1->p[n], U1_sum_size);
        if (iter > 0)
        {
            printf("\n\"LUU_lb\": %lf \t%ld   \t%ld    \t%ld   \t%ld", t_luu3,
                L3_nnz, L3_sum_size, U3_nnz, U3_sum_size);
        }
        printf("\nLUU: \t%lf \t%ld   \t%ld    \t%ld   \t%ld", t_luu,
            L_nnz, L2_sum_size, U_nnz, U2_sum_size);
        printf("\nLLU/LUU: %.4lf \t%.4lf \t%.4lf    \t%.4lf   \t%.4lf\n\n",
            t_llu/t_luu, (double)(L1->p[n])/(double)(L_nnz),
            (double)L1_sum_size/(double)L2_sum_size,
            (double)(U1->p[n])/(double)(U_nnz),
            (double)(U1_sum_size)/(double)(U2_sum_size));

        if (iter > 0)
        {
            printf("\nSPEX Left LU Factorization time: \t%lf", t1);
            printf("\n\"lower bound\" of SPEX LU Update time: \t%lf", t3);
            printf("\nSPEX LU Update time: \t\t\t%lf\n\n", t2);
            printf("\ntime of Left LU/ time of update: \t%lf\n\n", t1/t2);
        }

        iter++;
    }

    FREE_WORKSPACE;
    return 0;
}

