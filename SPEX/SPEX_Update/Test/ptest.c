//------------------------------------------------------------------------------
// SPEX_Update/Test/ptest.c: performance test for SPEX_Update
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

/*
 * performance test
 * ./ptest
 *           run ptest with lp problem lp_afiro, files located in
 *           TestMats/LPnetlib/lp_afiro/lp_afiro_*.mtx, no init basis
 *
 * ./ptest lp_25fv47
 *           run ptest with lp problem lp_25fv47, files located in
 *           TestMats/LPnetlib/lp_25fv47/lp_25fv47_*.mtx, no init basis
 *
 * ./ptest lp_25fv47 0
 *           run ptest with lp problem lp_25fv47, files located in
 *           TestMats/LPnetlib/lp_25fv47/lp_25fv47_*.mtx, no init basis
 *
 * ./ptest lp_25fv47 1 // does not work correctly
 *           run ptest with lp problem lp_25fv47, files located in
 *           TestMats/LPnetlib/lp_25fv47/lp_25fv47_*.mtx with init basis in
 *           the same directory
 */

// uncomment the following if basis info in each simplex iteration is wished
// to be written to an output file.
#define MY_PRINT_BASIS

#define FREE_WORKSPACE                           \
{                                                \
    glp_delete_prob(LP);                         \
    glp_free_env();                              \
    if (out_file != NULL) fclose(out_file);      \
    if (basis_file != NULL) fclose(basis_file);  \
    fclose(result_file);                         \
    SPEX_matrix_free(&Prob_A, option);           \
    SPEX_matrix_free(&Prob_c, option);           \
    SPEX_matrix_free(&tempA, option);            \
    SPEX_matrix_free(&A_CSC, option);            \
    SPEX_matrix_free (&x1, option);              \
    SPEX_matrix_free(&A_DCSC, option);           \
    SPEX_matrix_free(&b, option);                \
    SPEX_matrix_free(&b_dbl, option);            \
    SPEX_matrix_free(&c, option);                \
    SPEX_matrix_free(&c_new, option);            \
    SPEX_matrix_free(&basic_sol, option);        \
    SPEX_matrix_free(&y, option);                \
    SPEX_matrix_free(&y_sol, option);            \
    SPEX_matrix_free(&vk, option);               \
    SPEX_factorization_free(&F1, option);        \
    SPEX_factorization_free(&F2, option);        \
    SPEX_symbolic_analysis_free(&analysis, option);\
    SPEX_FREE(option);                           \
    mpq_clear(obj); mpz_clear(tmpz);             \
    mpq_clear(minq); mpq_clear(tmpq1);           \
    mpq_clear(tmpq2);mpq_clear(maxq);            \
    SPEX_FREE(basis);                            \
    SPEX_FREE(used_as_basis);                    \
    SPEX_FREE(col_val);                          \
    SPEX_FREE(col_ind);                          \
    SPEX_finalize() ;                            \
}

#include "test.h"
#include <assert.h>

int main( int argc, char* argv[])
{
    //char *prob_name = "lp_80bau3b";
    //char *prob_name = "lp_25fv47";//5.5018458883E+03
    char *prob_name = "lp_afiro"; // optimal: -4.6475314E+02
    //char *prob_name = "aa5";
    int init_basis = 0;
    if (argc >= 2)
    {
        prob_name = argv[1];
    }
    if (argc >= 3)
    {
        init_basis = atoi(argv[2]);
    }
    SPEX_info info;
    //------------------------------------------------------------------
    // Initialize for SPEX libray
    //------------------------------------------------------------------

    SPEX_initialize () ;

    //------------------------------------------------------------------
    // Allocate memory
    //------------------------------------------------------------------

    int64_t n, i, j, k, p, nz, nvars;
    int sgn;
    double t1 = 0, t2 = 0, t3 = 0, t_solve = 0;
    double z0 = 0;
    SPEX_options option = NULL;
    SPEX_matrix Prob_A = NULL, Prob_c = NULL, tempA = NULL, b_dbl = NULL;
    SPEX_matrix A_CSC = NULL, x1 = NULL, A_DCSC = NULL;
    SPEX_matrix b = NULL, c = NULL,  basic_sol = NULL, y = NULL,
                c_new = NULL, y_sol = NULL;
    SPEX_matrix vk = NULL;
    SPEX_factorization *F1 = NULL, *F2 = NULL;
    mpz_t tmpz;
    mpq_t obj, minq, maxq, tmpq1, tmpq2;
    SPEX_symbolic_analysis* analysis = NULL;
    int64_t *basis = NULL, *used_as_basis = NULL;
    double *col_val = NULL;
    int *col_ind = NULL;
    clock_t start_llu, start_luu, end_llu, end_luu, start1, end1, start2, end2,
            start_luu3 = 0, end_luu3 = 0;
    FILE *out_file = NULL;
    FILE *result_file = NULL;
    FILE *basis_file = NULL;
    char file_name[1000];
    glp_prob *LP;
    glp_smcp parm;
    LP = glp_create_prob();
    glp_init_smcp(&parm);
    parm.it_lim = 1;
    OK(SPEX_create_default_options(&option));
    OK(SPEX_mpz_init(tmpz));
    OK(SPEX_mpq_init(obj));
    OK(SPEX_mpq_init(minq));
    OK(SPEX_mpq_init(maxq));
    OK(SPEX_mpq_init(tmpq1));
    OK(SPEX_mpq_init(tmpq2));

    //--------------------------------------------------------------------------
    // open output file
    //--------------------------------------------------------------------------
#ifdef MY_PRINT_BASIS
    out_file = fopen("result.out", "w");
    if (out_file == NULL)
    {
        printf("Could not open file!\n");
        FREE_WORKSPACE;
        return 0;
    }
#endif

    sprintf(file_name, "Results/LPnetlib/%s_chol.txt",prob_name);
    result_file = fopen(file_name, "w");
    if (result_file == NULL)
    {
        printf("Could not open file: %s!\n", file_name);
        FREE_WORKSPACE;
        return 0;
    }

    //--------------------------------------------------------------------------
    // read in matrix
    //--------------------------------------------------------------------------
    sprintf(file_name, "TestMats/LPnetlib/%s/%s",prob_name,prob_name);
    // read matrix A, arrays b, lb, ub and c to LP
    z0 = 0;
    OK(SPEX_construct_LP(LP, &Prob_A, &b_dbl, &Prob_c, &z0, file_name, option));
    n = Prob_A->m;
    nvars = Prob_A->n;

    //--------------------------------------------------------------------------
    // generate initial basis matrix
    //--------------------------------------------------------------------------
    // allocate A_DCSC with n sparse vectors with initially 0 nnz
    OK(SPEX_matrix_allocate(&A_DCSC, SPEX_DYNAMIC_CSC, SPEX_MPZ, n, n, 0,
        false, true, option));
    OK(SPEX_mpq_set(A_DCSC->scale, Prob_A->scale))
    basis = (int64_t*) SPEX_malloc(n * sizeof(int64_t));
    used_as_basis = (int64_t*) SPEX_malloc(nvars * sizeof(int64_t));
    if (!basis || !used_as_basis)
    {
        FREE_WORKSPACE;
        return 0;
    }

    // if initial basis is available, read in via file
    if (init_basis == 1)
    {
        sprintf(file_name, "TestMats/LPnetlib/%s/%s_basis.txt",prob_name,
            prob_name);
        printf("reading basis from file %s\n",file_name);
        basis_file = fopen(file_name,"r");
        if (basis_file == NULL)
        {
            printf("Could not open file: %s!\n", file_name);
            FREE_WORKSPACE;
            return 0;
        }
        // get basis indices and build basis matrix
        for (i = 0; i < n; i++)
        {
            if (fscanf(basis_file, "%ld, ", &j) == EOF)
            {
                printf("Error while reading file: %s!\n", file_name);
                FREE_WORKSPACE;
                return 0;
            }
            basis[i] = j;

            if (j < n)
            {
                mpq_get_den(tmpz, Prob_A->scale);
                int r;
                OK(SPEX_mpz_cmp_ui(&r, tmpz, 1));
                assert(r == 0);// denominator must be 1
                if (A_DCSC->v[i]->nzmax < 1)
                {
                    OK(SPEX_vector_realloc(A_DCSC->v[i], 1, option));
                }
                mpq_get_num(A_DCSC->v[i]->x[0], Prob_A->scale);
                A_DCSC->v[i]->i[0] = j;
                A_DCSC->v[i]->nz = 1;
            }
            else
            {
                j = j-n;
                nz = Prob_A->p[j+1]-Prob_A->p[j];
                if (A_DCSC->v[i]->nzmax < nz)
                {
                    OK(SPEX_vector_realloc(A_DCSC->v[i], nz, option));
                }
                nz = 0;
                for (p = Prob_A->p[j]; p < Prob_A->p[j+1]; p++)
                {
                    A_DCSC->v[i]->i[nz] = Prob_A->i[p];
                    OK(SPEX_mpz_set(A_DCSC->v[i]->x[nz], Prob_A->x.mpz[p]));
                    nz++;
                }
                A_DCSC->v[i]->nz = nz;
            }
        }

        // get a CSC x MPZ A_CSC from dynamic_CSC x MPZ A_DCSC
        OK(SPEX_matrix_copy(&A_CSC, SPEX_CSC, SPEX_MPZ, A_DCSC, option));

        //------------------------------------------------------------------
        // perform LU factorization for the initial matrix A_CSC
        //------------------------------------------------------------------
        //start_llu = clock();

        // perform symbolic analysis by getting the column preordering of A
        OK(SPEX_lu_analyze(&analysis, A_CSC, option));

        // perform the SPEX Left LU factorization to obtain matrices L
        // and U and a row permutation P such that PAQ = LDU.
        OK(SPEX_lu_factorize(&F1, A_CSC, analysis, option));

        // create a mpz matrix of b
        OK(SPEX_matrix_copy(&b, SPEX_DENSE, SPEX_MPZ, b_dbl, option));
    }
    else
    {
        glp_adv_basis(LP, 0);
        //glp_std_basis(LP);
        //glp_cpx_basis(LP);
        glp_factorize(LP); // in order to use glp_get_bhead

        printf("getting initial basic variables....\n");
        parm.msg_lev = GLP_MSG_OFF;
        printf("OPT(%d) FEAS(%d) INFEAS(%d) NOFEAS(%d) UNBND(%d) UNDEF(%d)\n",
            GLP_OPT, GLP_FEAS, GLP_INFEAS, GLP_NOFEAS, GLP_UNBND, GLP_UNDEF);
        int64_t count=0;
        double last_obj = glp_get_obj_val(LP);
        while (glp_get_status(LP) != GLP_FEAS)
        {
            if (glp_get_obj_val(LP) >= last_obj)
            {
                parm.it_lim++;
            }
            last_obj = glp_get_obj_val(LP);
            glp_simplex(LP, &parm);
            printf("%ld: %d, obj: %lf\n",count++,
                glp_get_status(LP), glp_get_obj_val(LP)+z0);
        }
        parm.it_lim = 1;

        col_val =  (double*) SPEX_malloc((n+1)*sizeof(double));
        col_ind =  (int*)    SPEX_malloc((n+1)*sizeof(int));
        if (!col_val  || !col_ind)
        {
            FREE_WORKSPACE;
            return 0;
        }

        int num_of_new_basis = -1;
        last_obj = glp_get_obj_val(LP)+z0;

        while (glp_get_status(LP) != GLP_OPT)
        {
            // check if there is any non-basic variable that has nonzero value
            // printf("BS(%d) NLOWB(%d) NUPB(%d) NFREE(%d) NFIX(%d)\n",
            //    GLP_BS, GLP_NL, GLP_NU, GLP_NF, GLP_NS);
            int illegal_count = 0;
            for (j = 0; j < nvars; j++)
            {
                double lb = glp_get_col_lb(LP, j+1);
                double ub = glp_get_col_ub(LP, j+1);
                assert (lb == 0);
                int col_stat = glp_get_col_stat(LP, j+1);
                double xj   = glp_get_col_prim(LP, j+1);
                if (glp_get_col_bind(LP, j+1) == 0 && xj != 0)
                {
                    printf("%d: x[%ld]=%lf is nonbasic, bnd=[%lf %lf]"
                        " state=%d\n",
                        illegal_count++, j,xj,lb,ub, col_stat);
                    assert(col_stat == GLP_NU);

                    // update LP in glpk
                    int colj_nz = glp_get_mat_col(LP, j+1, col_ind, col_val);
                    for (k = 1; k <= colj_nz; k++)
                    {
                        i = col_ind[k]-1;
                        b_dbl->x.fp64[i] -= (double) (col_val[k]*ub);
                        glp_set_row_bnds(LP, i+1, GLP_FX, b_dbl->x.fp64[i],0.0);
                        col_val[k] *= -1;
                    }
                    glp_set_mat_col(LP, j+1, colj_nz, col_ind, col_val);
                    z0 += glp_get_obj_coef(LP,j+1)*ub;
                    glp_set_obj_coef(LP,j+1, -1*glp_get_obj_coef(LP,j+1));
                    glp_set_col_stat(LP, j+1, GLP_NL);

                    // update the SPEX matrices
                    for (p = Prob_A->p[j]; p < Prob_A->p[j+1]; p++)
                    {
                        OK(SPEX_mpz_neg(Prob_A->x.mpz[p], Prob_A->x.mpz[p]));
                    }
                    OK(SPEX_mpz_neg(Prob_c->x.mpz[j], Prob_c->x.mpz[j]));
                }
            }
            if (illegal_count > 0)
            {
                glp_simplex(LP, &parm);
                printf("%ld: %d, obj: %lf\n",count++,
                    glp_get_status(LP), glp_get_obj_val(LP));
                continue;
            }

            // create a mpz matrix of b
            OK(SPEX_matrix_free(&b, option));
            OK(SPEX_matrix_copy(&b, SPEX_DENSE, SPEX_MPZ, b_dbl, option));

            bool one_more_simplex = false;

            // free memory allocated before allocating new memory
            OK(SPEX_matrix_free(&A_CSC, option));
            OK(SPEX_factorization_free(&F1, option));
            OK(SPEX_matrix_free (&x1, option));
            SPEX_symbolic_analysis_free(&analysis, option);

            // get basis indices and build basis matrix
            for (i = 0; i < n; i++)
            {
                j = glp_get_bhead(LP, i+1)-1;
                if (num_of_new_basis == -1)
                {
                    basis[i] = j;
                }
                else
                {
                    if (basis[i] != j)
                    {
                        basis[i] = j;
                        num_of_new_basis++;
                    }
                }

                if (j < n)
                {
                    mpq_get_den(tmpz, Prob_A->scale);
                    int r;
                    OK(SPEX_mpz_cmp_ui(&r, tmpz, 1));
                    assert(r == 0);// denominator must be 1
                    if (A_DCSC->v[i]->nzmax < 1)
                    {
                        OK(SPEX_vector_realloc(A_DCSC->v[i], 1, option));
                    }
                    mpq_get_num(A_DCSC->v[i]->x[0], Prob_A->scale);
                    A_DCSC->v[i]->i[0] = j;
                    A_DCSC->v[i]->nz = 1;
                }
                else
                {
                    j = j-n;
                    if (glp_get_col_prim(LP, j+1) < -1e-6)
                    {
                        printf("x[%ld]=%f\n",j,glp_get_col_prim(LP,j+1));
                    }
                    nz = Prob_A->p[j+1]-Prob_A->p[j];
                    if (A_DCSC->v[i]->nzmax < nz)
                    {
                        OK(SPEX_vector_realloc(A_DCSC->v[i], nz, option));
                    }
                    nz = 0;
                    for (p = Prob_A->p[j]; p < Prob_A->p[j+1]; p++)
                    {
                        A_DCSC->v[i]->i[nz] = Prob_A->i[p];
                        OK(SPEX_mpz_set(A_DCSC->v[i]->x[nz], Prob_A->x.mpz[p]));
                        nz++;
                    }
                    A_DCSC->v[i]->nz = nz;
                }
            }

            // get a CSC x MPZ A_CSC from dynamic_CSC x MPZ A_DCSC
            OK(SPEX_matrix_copy(&A_CSC, SPEX_CSC, SPEX_MPZ, A_DCSC, option));

            //------------------------------------------------------------------
            // perform LU factorization for the initial matrix A_CSC
            //------------------------------------------------------------------
            //start_llu = clock();

            // perform symbolic analysis by getting the column preordering of A
            OK(SPEX_lu_analyze(&analysis, A_CSC, option));

            // perform the SPEX Left LU factorization to obtain matrices L
            // and U and a row permutation P such that PAQ = LDU.
            OK(SPEX_lu_factorize(&F1, A_CSC, analysis, option));

            //end_llu = clock();

            //------------------------------------------------------------------
            // Solve LDU x = b
            //------------------------------------------------------------------

            OK(SPEX_lu_solve(&x1, F1, b, option));

            // check if any solution is infeasible
            double max_diff = 0, diff;
            for (i = 0; i < n; i++)
            {
                j = basis[i];
                if (j >= n)
                {
                    j= j-n;
                    double xi    = glp_get_col_prim(LP, j+1);
                    diff = fabs(mpq_get_d(x1->x.mpq[i])-xi);
                    if (diff > 1e-5)
                    {
                        if (diff > max_diff) max_diff = diff;
                        printf("big difference in solution! x[%ld]=%lf!=%lf\n",
                            j, mpq_get_d(x1->x.mpq[i]), xi);
                    }
                    /*else
                    {
                        printf("x[%ld]=%lf==%lf\n",
                            j, mpq_get_d(x1->x.mpq[i]), xi);
                    }*/

                    if (mpq_get_d(x1->x.mpq[i]) < 0)
                    {
                        gmp_printf("infeasible solution: x[%ld]=%Qd\n",j,
                            x1->x.mpq[i]);
                        one_more_simplex = true;
                        break;
                    }
                }
            }
            if (max_diff > 0)
            {
                printf("max difference in solution is %f\n", max_diff);
            }

            // perform another iteration of simplex if needed
            if (one_more_simplex)
            {
                if (num_of_new_basis == 0 ||
                    (num_of_new_basis != -1 &&
                     glp_get_obj_val(LP)+z0 >= last_obj))
                {
                    parm.it_lim++;
                }
                last_obj = glp_get_obj_val(LP)+z0;
                glp_simplex(LP, &parm);
                printf("%ld: %d, obj: %lf, new basis count: %d\n",count++,
                    glp_get_status(LP), glp_get_obj_val(LP)+z0,
                    num_of_new_basis);
                num_of_new_basis = 0;
            }
            else
            {
                break;
            }
        }
        SPEX_FREE(col_val);
        SPEX_FREE(col_ind);
        printf("preprocess finished\n");

        if (glp_get_status(LP) == GLP_OPT)
        {
            printf("LP is optimized!\n");
            FREE_WORKSPACE;
            return 0;
        }
    }

#if 0
    for (i = 0; i < n; i++)
    {
        j = basis[i];
        if(j < n)
        {
            printf("%ld ", j);
            // These rows are independent from other variables, update
            // corresponding entries in b to make sure they get correct results.
            OK(SPEX_mpq_get_den(tmpz, Prob_A->scale));
            OK(SPEX_mpz_mul(b->x.mpz[j], b->x.mpz[j], tmpz));
        }
    }

    count=0;
    int num_basis=0;
    double glpk_obj=0.0,obj2=0.0,obj3=0.0,obj4=0.0;
    /*for(i = 0; i<n; i++)
    {

        if (glp_get_row_prim(LP,i+1)!=0)
        printf("%ld: %lf(%ld)\n",count++,glp_get_row_prim(LP,i+1),i);
    }*/
    for(i = 0; i<Prob_A->n; i++)
    {
        if (glp_get_col_prim(LP,i+1)!=0)
        {
             OK(SPEX_mpq_set_z(tmpq1, Prob_c->x.mpz[i]));
             OK(SPEX_mpq_div(tmpq1, tmpq1, Prob_c->scale));
             if(fabs(mpq_get_d(tmpq1)-glp_get_obj_coef(LP,i+1))>1e-5)
                 printf("%lf     ",mpq_get_d(tmpq1)-glp_get_obj_coef(LP,i+1));

    //         printf("%ld:%lf(%ld) coef:%lf==%lf %d\n",count++,
    //             glp_get_col_prim(LP,i+1),i,mpq_get_d(tmpq1),
    //             glp_get_obj_coef(LP,i+1),glp_get_col_bind(LP,i+1));
             glpk_obj+=glp_get_col_prim(LP,i+1)*glp_get_obj_coef(LP,i+1);
             obj2+=glp_get_col_prim(LP,i+1)*mpq_get_d(tmpq1);
             if (glp_get_col_bind(LP,i+1)!=0)
             {
                 obj3+=glp_get_col_prim(LP,i+1)*glp_get_obj_coef(LP,i+1);
             }
        }
         if (glp_get_col_bind(LP,i+1)!=0)
         {
             obj4+=glp_get_col_prim(LP,i+1)*glp_get_obj_coef(LP,i+1);
             num_basis++;
         }
    }
    printf("\ninit basic variabls generated, corresponding obj value is \n%lf\n",
        glp_get_obj_val(LP));
    printf("%lf \n%lf\n%lf\n%lf(%d)\n", glpk_obj,obj2,obj3,obj4,num_basis);
    printf("z0=%lf\n",z0);
#endif

    //--------------------------------------------------------------------------
    // generate initial inputs for LU update
    //--------------------------------------------------------------------------
    for (j = 0; j < nvars; j++)
    {
        used_as_basis[j] = -1;
    }
    for (i = 0; i < n; i++)
    {
        j = basis[i];
        if (j >= n)   used_as_basis[j-n] = i;
    }
    // allocate space for vk
    OK(SPEX_matrix_allocate(&vk, SPEX_DYNAMIC_CSC, SPEX_MPZ, n, 1, 0, false,
        true, option));

    // allocate space for c and y
    OK(SPEX_matrix_allocate(&c, SPEX_DENSE, SPEX_MPZ, n, 1, n, false, true,
        option));
    OK(SPEX_matrix_allocate(&y, SPEX_DENSE, SPEX_MPZ, n, 1, n, false, true,
        option));

    int64_t new_col = 0, iter = 0;
    k = 0;
    while (iter < 100)
    {
        int64_t L_nnz = 0, U_nnz = 0, L3_nnz = 0, U3_nnz = 0;
        size_t U1_sum_size = 0, U2_sum_size = 0, U3_sum_size = 0,
               L1_sum_size = 0, L2_sum_size = 0, L3_sum_size = 0, s;
        t_solve = 0;
        start1 = clock();
        start2 = clock();
        // solve for x_basic
        OK(SPEX_matrix_free (&basic_sol, option));
        OK(SPEX_update_solve(&basic_sol, F1, b, option));
        end2 = clock();
        t_solve += (double) (end2 - start2) / CLOCKS_PER_SEC;

        // reset objective value = 0
        OK(SPEX_mpq_set_ui(obj, 0, 1));
        bool checked = false;
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

                // compute objective value
                OK(SPEX_mpq_set_z(tmpq1, c->x.mpz[i]));
                OK(SPEX_mpq_mul(tmpq1, tmpq1, basic_sol->x.mpq[i]));
                OK(SPEX_mpq_add(obj, obj, tmpq1));

#if 1
                if (mpq_sgn(basic_sol->x.mpq[i]) < 0)
                {
                    printf("x[%ld]=%lf<0\n",j,mpq_get_d(basic_sol->x.mpq[i]));
                    //gmp_printf("exact x[%ld]=%Qd\n",j,basic_sol->x.mpq[i]);
                    //OK(SPEX_PANIC);
                    if (!checked)
                    {
                        checked = true;
                        bool Is_correct;
                        OK(MY_update_verify(&Is_correct, F1, A_CSC, option));
                        assert(Is_correct);
                    }
                }
#else
                SPEX_gmp_printf("%ld xz[%ld]= %Qd \n",i,j,basic_sol->x.mpq[i]);
#endif
            }
        }

        // write current basis to output file
#ifdef MY_PRINT_BASIS
        for (i = 0 ; i < n; i++)
        {
            fprintf(out_file, "%ld, ",basis[i]);
        }
        fprintf(out_file,"\n");
#endif

        // set c->scale = Prob_c->scale
        OK(SPEX_mpq_set(c->scale, Prob_c->scale));

        // compute the real objective value with c->scale applied
        OK(SPEX_mpq_div(obj, obj, Prob_c->scale));
        //OK(SPEX_gmp_printf("obj value = %Qd\n",obj));
        printf("obj value = %f\n", mpq_get_d(obj)+z0);

        //----------------------------------------------------------------------
        // find the entering variable
        //----------------------------------------------------------------------
        // solve A'*c_new = c for updated coefficient for objective function
        SPEX_matrix_free(&c_new, option);
        start2 = clock();
        OK(SPEX_update_tsolve(&c_new, F1, c, option));
        end2 = clock();
        t_solve += (double) (end2 - start2) / CLOCKS_PER_SEC;

        new_col = -1;
        // find the most positive coefficient since we want to min obj
        OK(SPEX_mpq_set_ui(maxq, 0, 1));// maxq = 0
        for (j = 0; j < Prob_A->n; j++)
        {
            if (used_as_basis[j] >= 0) {continue;} // skip basis variables

            // compute -c[j]+c_new'*A(:,j)
            OK(SPEX_mpq_set_ui(tmpq1, 0, 1));
            // iterate across nnz in column j of A to get tmpq1 = c_new'*A(:,j)
            for (p = Prob_A->p[j]; p < Prob_A->p[j+1]; p++)
            {
                i = Prob_A->i[p];
                OK(SPEX_mpq_set_z(tmpq2, Prob_A->x.mpz[p]));
                OK(SPEX_mpq_mul(tmpq2, tmpq2, c_new->x.mpq[i]));
                OK(SPEX_mpq_add(tmpq1, tmpq1, tmpq2));
            }
            OK(SPEX_mpq_div(tmpq1, tmpq1, Prob_A->scale));
            // tmpq2 = c[j]
            OK(SPEX_mpq_set_z(tmpq2, Prob_c->x.mpz[j]));
            OK(SPEX_mpq_div(tmpq2, tmpq2, Prob_c->scale));
            // tmpq1 = tmpq1 - tmpq2
            mpq_sub(tmpq1, tmpq1, tmpq2);
            // use exact comparison
            OK(SPEX_mpq_cmp(&sgn, tmpq1, maxq));
            if (sgn > 0)
            {
                OK(SPEX_mpq_set(maxq, tmpq1));
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
        vk->v[0]->nz = Prob_A->p[new_col+1]-Prob_A->p[new_col];
        if (vk->v[0]->nzmax < vk->v[0]->nz)
        {
            OK(SPEX_vector_realloc(vk->v[0], vk->v[0]->nz, option));
        }
        i = 0;
        j = 0;
        // iterate across nnz in the entering column of A
        for (p = Prob_A->p[new_col]; p < Prob_A->p[new_col+1]; p++)
        {
            j = Prob_A->i[p];
            // dense y
            OK(SPEX_mpz_set(y->x.mpz[j], Prob_A->x.mpz[p]));

            // sparse vk
            vk->v[0]->i[i] = j;
            OK(SPEX_mpz_set(vk->v[0]->x[i], Prob_A->x.mpz[p]));
            i++;
        }
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

        // solve for Ay_sol = y
        SPEX_matrix_free(&y_sol, option);
        start2 = clock();
        OK(SPEX_update_solve(&y_sol, F1, y, option));
        end2 = clock();
        t_solve += (double) (end2 - start2) / CLOCKS_PER_SEC;

        // perform ratio test to find existing variable
        k = -1;
        for (i = 0; i < n; i++)
        {
            if (mpq_sgn(y_sol->x.mpq[i]) > 0)
            {
                if (mpq_sgn(basic_sol->x.mpq[i]) == 0)
                {
                    k = i;
                    printf("found 0!! This update won't improve obj value\n");
                    break;
                }

                // basic_sol[i] = basic_sol[i]/y_sol[i] (it's ok to modify
                // basic_sol here since it will be deleted/updated in next loop)
                OK(SPEX_mpq_div(basic_sol->x.mpq[i], basic_sol->x.mpq[i],
                                y_sol->x.mpq[i]));
                // printf(" x/y=%f \n",mpq_get_d(basic_sol->x.mpq[i]));
                if (k == -1)
                {
                    OK(SPEX_mpq_set(minq, basic_sol->x.mpq[i]));
                    k = i;
                }
                else
                {
                    OK(SPEX_mpq_cmp(&sgn, basic_sol->x.mpq[i], minq));
                    if (sgn < 0)
                    {
                        /*
                        OK(SPEX_gmp_printf("%Qd <%Qd\n",basic_sol->x.mpq[i],
                                           minq));
                        printf("%f < %f\n",mpq_get_d(basic_sol->x.mpq[i]),
                                           mpq_get_d(minq));
                        */
                        OK(SPEX_mpq_set(minq, basic_sol->x.mpq[i]));
                        k = i;
                    }
                }
            }
        }
        if (k == -1)
        {
            printf("LP is unbound!\n");
            FREE_WORKSPACE;
            return 0;
        }
        end1 = clock();

        // check basis
        /*
        for (i = 0; i < n; i++)
        {
            if(used_as_basis[i] == -1) continue;
            printf("\n%ld:\n",i);
            for (j = i+1; j < n; j++)
            {
                printf ("[%ld %ld]",used_as_basis[i], used_as_basis[j]);
                if (used_as_basis[i] == used_as_basis[j])
                {
                    printf("%ld %ld->%ld, basis[%ld]=%ld\n",i,j,
                         used_as_basis[i], used_as_basis[i],
                         basis[used_as_basis[i]]);
                    abort();
                }
            }
        }*/

        // print results
        printf("\nSolving 3 linear equations time: \t%lf", t_solve);
        printf("\nSearch k and new_col time: \t\t%lf\n",
           (double) (end1 - start1) / CLOCKS_PER_SEC);
        fprintf(result_file, "%lf \t%lf ", t_solve,
           (double) (end1 - start1) / CLOCKS_PER_SEC);
#ifdef MY_PRINT_BASIS
        fprintf(out_file,"prev basis k(%ld): %ld; new basis k: %ld\n",
            k, basis[k]-n, new_col);
        fprintf(out_file,
            "\n---------------------------------------------------------\n");
        fprintf(out_file,
            "----------%ld: replacing k(%ld) with new_col(%ld)-----------\n",
            iter, k, new_col);
        fprintf(out_file,
            "-----------------------------------------------------------\n");
#endif
        printf("prev basis k(%ld): %ld; new basis k: %ld\n",
            k, basis[k]-n, new_col);
        used_as_basis[new_col] = k;
        if (basis[k] >= n)
        {
            assert(used_as_basis[basis[k]-n] == k);
            used_as_basis[basis[k]-n] = -1;
        }
        basis[k] = new_col+n;
        printf("\n---------------------------------------------------------\n");
        printf("----------%ld: replacing k(%ld) with new_col(%ld)-----------\n",
            iter, k, new_col);
        printf("-----------------------------------------------------------\n");

        //----------------------------------------------------------------------
        // perform continuous LU update
        //----------------------------------------------------------------------
        start_luu = clock();

        OK(SPEX_update_lu_colrep(F1, vk, k, option));

        end_luu = clock();

        //----------------------------------------------------------------------
        // perform LU update using the direct factorization from last loop
        //----------------------------------------------------------------------

        if (iter > 0)
        {
            // perform conversion first to exclude the time for this step
            OK(SPEX_factorization_convert(F2, true, option));

            start_luu3 = clock();

            OK(SPEX_update_lu_colrep(F2, vk, k, option));

            end_luu3 = clock();

            for (i = 0; i < n; i++)
            {
                L3_nnz += F2->L->v[i]->nz;
                for (p = 0; p < F2->L->v[i]->nz; p++)
                {
                    OK(SPEX_mpz_sizeinbase(&s, F2->L->v[i]->x[p], 2));
                    L3_sum_size += s;
                }
                U3_nnz += F2->U->v[i]->nz;
                for (p = 0; p < F2->U->v[i]->nz; p++)
                {
                    OK(SPEX_mpz_sizeinbase(&s, F2->U->v[i]->x[p], 2));
                    U3_sum_size += s;
                }
            }
        }

        //----------------------------------------------------------------------
        // generate new matrix with vk inserted
        //----------------------------------------------------------------------
        OK(SPEX_update_matrix_colrep(A_DCSC, vk, k, option));
        OK(SPEX_matrix_free(&A_CSC, option));
        OK(SPEX_matrix_copy(&A_CSC, SPEX_CSC, SPEX_MPZ, A_DCSC, option));

        //----------------------------------------------------------------------
        // perform direct LU factorization for matrix A_CSC
        //----------------------------------------------------------------------
        SPEX_factorization_free(&F2, NULL);
        SPEX_symbolic_analysis_free(&analysis, option);
        start_llu = clock();

        // perform symbolic analysis by getting the column preordering of A
        OK(SPEX_lu_analyze(&analysis, A_CSC, option));

        // Now we perform the SPEX Left LU factorization to obtain matrices L
        // and U and a row permutation P such that PAQ = LDU.
        OK(SPEX_lu_factorize(&F2, A_CSC, analysis, option));
//        if (info == SPEX_OK) {printf("matrix is not singular!\n");}

        end_llu = clock();



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
            L_nnz += F1->L->v[i]->nz;
            for (p = 0; p < F1->L->v[i]->nz; p++)
            {
                OK(SPEX_mpz_sizeinbase(&s, F1->L->v[i]->x[p], 2));
                L2_sum_size += s;
            }
            U_nnz += F1->U->v[i]->nz;
            for (p = 0; p < F1->U->v[i]->nz; p++)
            {
                OK(SPEX_mpz_sizeinbase(&s, F1->U->v[i]->x[p], 2));
                U2_sum_size += s;
            }
        }
        for (p = 0; p < F2->L->p[n]; p++)
        {
            OK(SPEX_mpz_sizeinbase(&s, F2->L->x.mpz[p], 2));
            L1_sum_size += s;
        }
        for (p = 0; p < F2->U->p[n]; p++)
        {
            OK(SPEX_mpz_sizeinbase(&s, F2->U->x.mpz[p], 2));
            U1_sum_size += s;
        }
        printf("\n\ttime     \tL_nnz \tsum(Lb) \tU_nnz \t\tsum(Ub)\n");
        printf("\nLLU: \t%lf \t%ld   \t%ld    \t%ld   \t%ld", t_llu,
            F2->L->p[n], L1_sum_size, F2->U->p[n], U1_sum_size);
        fprintf(result_file, "\t%lf \t%ld \t%ld \t%ld \t%ld ", t_llu,
            F2->L->p[n], L1_sum_size, F2->U->p[n], U1_sum_size);
        if (iter > 0)
        {
            printf("\n\"LUU_lb\": %lf \t%ld   \t%ld    \t%ld   \t%ld", t_luu3,
                L3_nnz, L3_sum_size, U3_nnz, U3_sum_size);
            fprintf(result_file,"\t%lf \t%ld \t%ld \t%ld \t%ld ", t_luu3,
                L3_nnz, L3_sum_size, U3_nnz, U3_sum_size);
        }
        else
        {
            fprintf(result_file,"\t0.0 \t0 \t0 \t0 \t0 ");
        }
        printf("\nLUU: \t%lf \t%ld   \t%ld    \t%ld   \t%ld", t_luu,
            L_nnz, L2_sum_size, U_nnz, U2_sum_size);
        fprintf(result_file,"\t%lf \t%ld \t%ld \t%ld \t%ld\n", t_luu,
            L_nnz, L2_sum_size, U_nnz, U2_sum_size);
        printf("\nLLU/LUU: %.4lf \t%.4lf \t%.4lf    \t%.4lf   \t%.4lf\n\n",
            t_llu/t_luu, (double)(F2->L->p[n])/(double)(L_nnz),
            (double)L1_sum_size/(double)L2_sum_size,
            (double)(F2->U->p[n])/(double)(U_nnz),
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

