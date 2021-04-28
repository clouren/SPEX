//------------------------------------------------------------------------------
// SPEX_LU_Update/Test/ptest.c: performance test for SPEX_LU_Update library
//------------------------------------------------------------------------------

// SPEX_LU_Update: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_LU_Update/License for the license.

//------------------------------------------------------------------------------

/*
 * performance test
 */

#define FREE_WORKSPACE                           \
{                                                \
    if (mat_file!=NULL) fclose(mat_file);        \
    SPEX_FREE(option);                           \
    SPEX_finalize() ;                            \
}

#include "test.h"
//#include <time.h>
#include <assert.h>

int main( int argc, char* argv[])
{
    //char *prob_name = "lp_80bau3b";
    char *prob_name = "lp_25fv47";
    //char *prob_name = "aa5";
    if (argc >= 2)
    {
        prob_name = argv[1];
    }
    SPEX_info info;
    //------------------------------------------------------------------
    // Initialize SPEX CHOLMOD process
    //------------------------------------------------------------------

    SPEX_initialize () ;

    //------------------------------------------------------------------
    // Allocate memory
    //------------------------------------------------------------------

    int64_t n, i, j, p, nz;
    int sgn;
    double t1 = 0, t2 = 0;
    SPEX_options* option = NULL;
    SPEX_matrix *Prob_A = NULL, *Prob_b = NULL, *Prob_c;
    SPEX_matrix *L1 = NULL, *U1 = NULL, *rhos = NULL, *A1 = NULL;
    SPEX_mat *L2 = NULL, *U2 = NULL, *A2 = NULL, *b = NULL, *c = NULL,
             *x  = NULL, *y = NULL;
    SPEX_vector *tmpv, *vk = NULL;
    mpz_t *d = NULL, *sd = NULL;
    mpz_t z, minz, tmpz;
    mpq_t minq, tmpq, one;
    SPEX_matrix *S = NULL;
    SPEX_LU_analysis* analysis = NULL;
    int64_t *P1_inv = NULL, *P = NULL, *P_inv = NULL, *Q = NULL, *Q_inv = NULL;
    int64_t *basis = NULL, *h = NULL, *used_as_basis = NULL;
    clock_t start_llu, start_luu, end_llu, end_luu;
    FILE *mat_file = NULL;
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
    // read in matrix
    //--------------------------------------------------------------------------
    char file_name[1000] = "TestMats/";
    strcat(file_name, prob_name);
    // read matrix A, arrays b, lb, ub and c to LP
    SPEX_construct_LP(LP, &Prob_A, &Prob_b, &Prob_c, file_name, option);
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
    printf("%d %d %d %d\n", GLP_UNDEF, GLP_FEAS,GLP_INFEAS,GLP_NOFEAS);
    printf("%d %d %d %d %d %d\n", GLP_OPT, GLP_FEAS,GLP_INFEAS,GLP_NOFEAS,GLP_UNBND,GLP_UNDEF);
#if 1
    while (glp_get_status(LP) != GLP_FEAS)
    {
        glp_simplex(LP, &parm);
    }
    for (i = 0; i < n; i++)
    {
        basis[i] = glp_get_bhead(LP, i+1)-1;
        if (basis[i] < n && glp_get_mat_row(LP, basis[i]+1, NULL, NULL) != 0)
        {
            //printf("basis[%ld]=%ld\n",i,basis[i]);
            GOTCHA;
            return 0;;
        }
    }
    double opt_val =0.0;
    for (j = 0; j<Prob_A->n;j++)
    {
        if (glp_get_col_stat(LP,j+1)==GLP_BS)
        opt_val += glp_get_col_prim(LP,j+1)*glp_get_obj_coef(LP,j+1);
        if (glp_get_col_prim(LP,j+1)<0)
            printf("panic %ld\n",j);
    }
    printf("%f\n",opt_val);
#else
    while (glp_get_status(LP) != GLP_OPT)
    {
        // run one iteration of simplex and find the new basis
        glp_simplex(LP, &parm);
        bool found = false;
        for (i = 0; i < n; i++)
        {
            basis[i] = glp_get_bhead(LP, i+1)-1;
            if (basis[i] < n && glp_get_mat_row(LP, basis[i]+1, NULL, NULL) != 0)
            {
                //printf("basis[%ld]=%ld\n",i,basis[i]);
                found = true;
                break;
            }
        }
        //printf("%s\n",glp_get_status(LP)==GLP_OPT?"optimal":"not yet");
        //if (!found) {break;}
        printf("%d %d %d\n",glp_get_prim_stat(LP),glp_get_dual_stat(LP),glp_get_status(LP));
        double opt_val =0.0;
        for (j = 0; j<Prob_A->n;j++)
        {
            if (glp_get_col_stat(LP,j+1)==GLP_BS)
            opt_val += glp_get_col_prim(LP,j+1)*glp_get_obj_coef(LP,j+1);
            if (glp_get_col_prim(LP,j+1)<0 && !found)
                printf("panic\n");
        }
        printf("%f\n",opt_val);
    }
#endif
    if (glp_get_status(LP) == GLP_OPT)
    {
        printf("cannot find basic variable set without row variables!\n");
        return 0;
    }
    GOTCHA;

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
    OK(SPEX_mat_alloc(&A2, n, n, true));
    for (i = 0; i < n; i++)
    {
        if (basis[i] < n)
        {
            OK(SPEX_vector_realloc(A2->v[i], 1));
            mpq_get_num(A2->v[i]->x[0], Prob_A->scale);
            A2->v[i]->i[0] = basis[i];
            A2->v[i]->nz = 1;
        }
        else
        {
            j = basis[i]-n;
            used_as_basis[j] = i;
            OK(SPEX_vector_realloc(A2->v[i], Prob_A->p[j+1]-Prob_A->p[j]));
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
    OK(SPEX_mat_to_CSC(&A1, A2, NULL, true, option));
    GOTCHA;

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
    // generate permutation vectors P, Q, P_inv, Q_inv and vectors d, sd,
    // and create S
    Q     = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    P     = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    Q_inv = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    P_inv = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    d  = SPEX_create_mpz_array(n);
    sd = SPEX_create_mpz_array(n);
    h  = (int64_t*) SPEX_calloc(n, sizeof(int64_t));
    if (!d || !sd || !P || !Q || !P_inv || !Q_inv || !h)
    {
        FREE_WORKSPACE;
        return 0;
    }
    OK(SPEX_matrix_allocate(&S, SPEX_DENSE, SPEX_MPQ, 2, n, n,
        false, true, option));
    for (i = 0; i < n; i++)
    {
        Q[i] = analysis->q[i];
        P_inv[i] = P1_inv[i];
        P[P_inv[i]] = i;
        Q_inv[Q[i]] = i;
        OK(SPEX_mpz_set(sd[i], SPEX_1D(rhos, i, mpz)));
        OK(SPEX_mpz_set(d[i], sd[i]));
        OK(SPEX_mpq_set_ui(SPEX_2D(S, 0, i, mpq), 1, 1));
        OK(SPEX_mpq_set_ui(SPEX_2D(S, 1, i, mpq), 1, 1));
    }

    // convert the factorization to SPEX_mat to be used in the update process
    OK(SPEX_CSC_to_mat(&L2, P, true,  L1, option));
    OK(SPEX_mat_canonicalize(L2, P));
    OK(SPEX_CSC_to_mat(&U2, Q, false, U1, option));// U2 stored in row-wise
    OK(SPEX_mat_canonicalize(U2, Q));

    // allocate space for vk
    OK(SPEX_vector_alloc(&vk, 0, true));

    // get SPEX_mat b from SPEX_matrix Prob_b
    OK(SPEX_mat_alloc(&b, 1, n, true));
    SPEX_FREE(b->v[0]->i); // no need for dense vector
    b->v[0]->x = Prob_b->x.mpz;
    Prob_b->x.mpz = NULL;
    OK(SPEX_matrix_free(&Prob_b, option));

    // allocate space for c
    OK(SPEX_mat_alloc(&c, 1, n, true));
    SPEX_FREE(c->v[0]->i); // no need for dense vector
    c->v[0]->x = SPEX_create_mpz_array(n);
    // allocate space for y
    OK(SPEX_mat_alloc(&y, 1, n, true));
    SPEX_FREE(y->v[0]->i); // no need for dense vector
    y->v[0]->x = SPEX_create_mpz_array(n);

    int64_t k = 0, new_col = 0;
    while (1)
    {
    GOTCHA;
        // solve for x_basic
        OK(SPEX_solve(&x, b, L2, U2, A2->scale, false, h, S, (const mpz_t*)sd, P,
            Q_inv, option));
        printf("scale sgn=%d\n",mpq_sgn(x->scale));

        // reset objective value z = 0
        OK(SPEX_mpz_set_ui(z, 0));
        for (i = 0; i < n; i++)
        {
            j = basis[i];
            if (j < n)
            {
                //j = j-1;
                // build vector c
                OK(SPEX_mpz_set_ui(c->v[0]->x[i], 0));
            }
            else
            {
                j = j-n;
                // build vector c
                OK(SPEX_mpz_set(c->v[0]->x[i], Prob_c->x.mpz[j]));

                // compute objective value z
                OK(SPEX_mpz_addmul(z, c->v[0]->x[i], x->v[0]->x[i]));

                OK(SPEX_mpq_set_z(tmpq, x->v[0]->x[i]));
                OK(SPEX_mpq_div(tmpq, tmpq, x->scale));
                if (mpz_sgn(x->v[0]->x[i]) * mpq_sgn(x->scale) < 0)
                {
                    printf("xz[%ld]<0 %d %d %f \n",j,mpz_sgn(x->v[0]->x[i]),mpq_sgn(x->scale),mpq_get_d(tmpq));
                }
                if (i == k)
                {
                    printf("x[%ld] = %f\n",j,mpq_get_d(tmpq));
                }
#if 0
                if (glp_get_col_prim(LP,j+1) < 0)
                    printf("     x[%ld]<0 %f %f\n",j,glp_get_col_prim(LP,j+1),mpq_get_d(tmpq) );
#endif
                /*if (glp_get_col_prim(LP,j+1) - mpq_get_d(tmpq) > 0.001 ||
                    mpq_get_d(tmpq) < 0)
                    printf("%f != %f\n",glp_get_col_prim(LP,j+1),mpq_get_d(tmpq));*/
            }
        }
        OK(SPEX_mpq_set_z(tmpq, z));
        OK(SPEX_mpq_div(tmpq, tmpq, x->scale));
        OK(SPEX_mpq_div(tmpq, tmpq, Prob_c->scale));
        //OK(SPEX_gmp_printf("obj value = %Qd\n",tmpq));
        printf("obj value = %f\n", mpq_get_d(tmpq));

        // solve for updated coefficient for objective function
        OK(SPEX_solve(&c, c, U2, L2, A2->scale, true, h, S, (const mpz_t*) sd, Q,
            P_inv, option));

        // find the entering variable
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
                OK(SPEX_mpz_addmul(tmpz, Prob_A->x.mpz[p], c->v[0]->x[i]));
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
            if (tmpd < mind)
            {
                //OK(SPEX_gmp_printf("c:%Zd <%Zd\n",tmpz,minz));
                //OK(SPEX_mpz_set(minz, tmpz));
                printf("c: %f < %f new_col = %ld --> %ld\n",tmpd, mind,new_col,j);
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
            OK(SPEX_mpz_set_ui(y->v[0]->x[i], 0));
        }
        // allocate space for vk if needed
        vk->nz = Prob_A->p[new_col+1]-Prob_A->p[new_col];
        if (vk->nzmax < vk->nz)
        {
            OK(SPEX_vector_realloc(vk, vk->nz));
        }
        i = 0;
        j = 0;
        // iterate across nnz in the entering column of A
        for (p = Prob_A->p[new_col]; p < Prob_A->p[new_col+1]; p++)
        {
            for (; j < Prob_A->i[p]; j++)
            {
                OK(SPEX_mpz_set_ui(y->v[0]->x[j], 0));
            }
            OK(SPEX_mpz_set(y->v[0]->x[j], SPEX_1D(Prob_A, p, mpz)));// dense y

            vk->i[i] = j;
            OK(SPEX_mpz_set(vk->x[i], SPEX_1D(Prob_A, p, mpz)));// sparse vk
            i++;
            j++;

        }
        vk->nz = i;

        // solve for Ay = y
        OK(SPEX_solve(&y, y, L2, U2, A2->scale,false, h, S, (const mpz_t*) sd, P,
            Q_inv, option));
        // perform ratio test to find existing variable
        k = -1;
        for (i = 0; i < n; i++)
        {
            int x_sgn = mpz_sgn(x->v[0]->x[i]) * mpq_sgn(x->scale);
            //if (x_sgn < 0) OK(SPEX_PANIC);

            if (x_sgn > 0 && mpz_sgn(y->v[0]->x[i]) > 0)
            {
                OK(SPEX_mpq_set_num(tmpq, x->v[0]->x[i]));
                OK(SPEX_mpq_set_den(tmpq, y->v[0]->x[i]));
                OK(SPEX_mpq_canonicalize(tmpq));
                if (mpq_sgn(x->scale) < 0)
                {
                    OK(SPEX_mpq_neg(tmpq, tmpq));
                }
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
                        printf("%f <%f\n",mpq_get_d(tmpq),mpq_get_d(minq));
                        OK(SPEX_mpq_set(minq, tmpq));
                        k = i;
                    }
                }
            }
        }
        if (k == -1)
        {
            printf("unknown error!\n");
            return 0;
        }
        used_as_basis[new_col] = k;
        basis[k] = new_col+n;
        used_as_basis[k] = -1;
        printf("------------------------------------------------------------\n");
        printf("--------------replacing k(%ld) with new_col(%ld)------------\n",
            k, new_col);
        printf("------------------------------------------------------------\n");

        //----------------------------------------------------------------------
        // generate new matrix with vk inserted
        //----------------------------------------------------------------------
        tmpv = A2->v[k]; A2->v[k] = vk; vk = tmpv;
        OK(SPEX_matrix_free(&A1, option));
        OK(SPEX_mat_to_CSC(&A1, A2, NULL, true, option));
        tmpv = A2->v[k]; A2->v[k] = vk; vk = tmpv;

        //----------------------------------------------------------------------
        // perform LU factorization for matrix A1
        //----------------------------------------------------------------------
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
        // perform LU update for matrix A2->A1
        //----------------------------------------------------------------------
        /*OK(SPEX_gmp_printf("%Qd\n",SPEX_LUU_2D(S,1,k)));
        OK(SPEX_gmp_printf("%Qd\n",SPEX_LUU_2D(S,2,k)));
        OK(SPEX_gmp_printf("%Qd\n",SPEX_LUU_2D(S,3,k)));
        mpq_t tmpq; mpq_init(tmpq);
        mpz_t tmpz; mpz_init(tmpz);
        mpz_t mysd; mpz_init(mysd);
        mpq_mul(tmp_q,SPEX_LUU_2D(S,1,k),SPEX_LUU_2D(S,1,k));
        mpq_get_den(tmpz,tmpq);
        mpz_mul(mysd,L2->v[k]->x[P[k]],*/
        start_luu = clock();

        OK(SPEX_LUU(A2, L2, U2, sd, S, P, P_inv, Q, Q_inv, &vk, k, option));

        end_luu = clock();
        for(int ii =0; ii<n;ii++)
        {
            for (int jj =ii+1;jj<n;jj++)
            {
                if (Q[jj] == Q[ii])
                {
                    printf("Q %d %d\n",ii,jj);
                    return 0;
                }
                if (P[jj] == P[ii])
                {
                    printf("P %d %d\n",ii,jj);
                    return 0;
                }
            }
        }

        //----------------------------------------------------------------------
        // print results
        //----------------------------------------------------------------------
        // Timing stats
        double t_llu = (double) (end_llu-start_llu)/CLOCKS_PER_SEC;
        double t_luu = (double) (end_luu - start_luu) / CLOCKS_PER_SEC;
        t1 += t_llu;
        t2 += t_luu;
        printf("\nSPEX Left LU Factorization time: \t%lf", t_llu);
        printf("\nSPEX LU Update time: \t\t\t%lf\n\n", t_luu);
        printf("\nSPEX Left LU Factorization time: \t%lf", t1);
        printf("\nSPEX LU Update time: \t\t\t%lf\n\n", t2);
        printf("\ntime of Left LU/ time of update: \t%lf\n\n", t1/t2);
    }

    FREE_WORKSPACE;
    return 0;
}

