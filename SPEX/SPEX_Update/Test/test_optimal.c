//------------------------------------------------------------------------------
// SPEX_Update/Test/test_optimal.c: test if the optimal value is valid
//------------------------------------------------------------------------------

// SPEX_Update: (c) 2020-2023, Chris Lourenco, Jinhao Chen,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "test.h"

int main (int argc, char *argv[])
{
    //char *prob_name = "lp_80bau3b";
    //char *prob_name = "lp_25fv47";//5.5018458883E+03
    char *prob_name = "lp_afiro"; // optimal: -4.6475314E+02
    //char *prob_name = "aa5";
    if (argc >= 2)
    {
        prob_name = argv[1];
    }
    SPEX_info info;

    //------------------------------------------------------------------
    // Initialize for SPEX libray
    //------------------------------------------------------------------

    SPEX_initialize ();

    //------------------------------------------------------------------
    // Allocate memory
    //------------------------------------------------------------------

    int64_t n, i, j, k, p, nz, nvars;
    double z0 = 0;
    SPEX_options option = NULL;
    SPEX_matrix Prob_A = NULL, Prob_c = NULL, b_dbl = NULL;
    SPEX_matrix A_CSC = NULL, x1 = NULL, A_DCSC = NULL;
    SPEX_matrix b = NULL;
    SPEX_factorization F1 = NULL;
    mpz_t tmpz;
    mpq_t obj, minq, maxq, tmpq1, tmpq2;
    SPEX_symbolic_analysis analysis = NULL;
    int64_t *basis = NULL, *used_as_basis = NULL;
    double *col_val = NULL;
    int *col_ind = NULL;
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
    for (j = 0; j < nvars; j++)
    {
        used_as_basis[j] = -1;
    };

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
    while (glp_get_status(LP) != GLP_OPT)
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
    printf("LP is optimized!\n");

    col_val =  (double*) SPEX_malloc((n+1)*sizeof(double));
    col_ind =  (int*)    SPEX_malloc((n+1)*sizeof(int));
    if (!col_val  || !col_ind)
    {
        FREE_WORKSPACE;
        return 0;
    }

    // check if there is any non-basic variable that has nonzero value
    // printf("BS(%d) NLOWB(%d) NUPB(%d) NFREE(%d) NFIX(%d)\n",
    //    GLP_BS, GLP_NL, GLP_NU, GLP_NF, GLP_NS);
    int illegal_count = 0;
    for (j = 0; j < nvars; j++)
    {
        double lb = glp_get_col_lb(LP, j+1);
        double ub = glp_get_col_ub(LP, j+1);
        ASSERT (lb == 0);
        int col_stat = glp_get_col_stat(LP, j+1);
        double xj   = glp_get_col_prim(LP, j+1);
        if (glp_get_col_bind(LP, j+1) == 0 && xj != 0)
        {
            printf("%d: x[%ld]=%lf is nonbasic, bnd=[%lf %lf]"
                " state=%d\n",
                illegal_count++, j,xj,lb,ub, col_stat);
            ASSERT (col_stat == GLP_NU);

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
    }

    // create a mpz matrix of b
    OK(SPEX_matrix_free(&b, option));
    OK(SPEX_matrix_copy(&b, SPEX_DENSE, SPEX_MPZ, b_dbl, option));

    // free memory allocated before allocating new memory
    OK(SPEX_matrix_free(&A_CSC, option));
    OK(SPEX_factorization_free(&F1, option));
    OK(SPEX_matrix_free (&x1, option));
    SPEX_symbolic_analysis_free(&analysis, option);

    // get basis indices and build basis matrix
    for (i = 0; i < n; i++)
    {
        j = glp_get_bhead(LP, i+1)-1;
        basis[i] = j;
        printf("%ld, ", j);

        if (j < n)
        {
            mpq_get_den(tmpz, Prob_A->scale);
            int r;
            OK(SPEX_mpz_cmp_ui(&r, tmpz, 1));
            ASSERT (r == 0);// denominator must be 1
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
            used_as_basis[j] = i;
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
    printf("\n");

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
            }
            // compute objective value
            OK(SPEX_mpq_set_z(tmpq1, Prob_c->x.mpz[j]));
            OK(SPEX_mpq_mul(tmpq1, tmpq1, x1->x.mpq[i]));
            OK(SPEX_mpq_add(obj, obj, tmpq1));
        }
    }

    // compute the real objective value with c->scale applied
    OK(SPEX_mpq_div(obj, obj, Prob_c->scale));
    printf("exact obj value = %f\n", mpq_get_d(obj)+z0);
    if (max_diff > 0)
    {
        printf("max difference in solution is %f\n", max_diff);
    }


    FREE_WORKSPACE;
    return 0;
}
