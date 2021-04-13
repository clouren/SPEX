//------------------------------------------------------------------------------
// SPEX_LU_Update/Test/test.c: performance test for SPEX_LU_Update library
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
    SPEX_options* option = NULL;
    SPEX_matrix *Prob = NULL, *Prob_trip = NULL;
    SPEX_matrix *L1 = NULL, *U1 = NULL, *rhos = NULL, *A1 = NULL;
    SPEX_mat *L2 = NULL, *U2 = NULL, *A2 = NULL;
    SPEX_vector *tmpv, *vk = NULL;
    mpz_t *d = NULL, *sd = NULL;
    mpq_t *S = NULL;
    SPEX_LU_analysis* analysis = NULL;
    int64_t *P1_inv = NULL;
    int64_t *P = NULL, *P_inv = NULL, *Q = NULL, *Q_inv = NULL, *basis = NULL;
    clock_t start_llu, start_luu, end_llu, end_luu;
    FILE *mat_file = NULL;
    glp_prob *LP;
    glp_smcp parm;
    LP = glp_create_prob();
    glp_init_smcp(&parm);
    parm.it_lim = 1;

    //--------------------------------------------------------------------------
    // read in matrix
    //--------------------------------------------------------------------------
    char file_name[1000] = "TestMats/";
    strcat(file_name, prob_name);
    char *suffix = ".mtx";
    strcat(file_name, suffix);
    mat_file = fopen(file_name, "r");
    if (mat_file == NULL)
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }
    OK(SPEX_create_default_options(&option));

    // Read matrix from given file as a triplet matrix
    OK(SPEX_mmread(&Prob_trip, mat_file, option));
    fclose(mat_file); mat_file = NULL;
    // load matrix to LP
    OK(SPEX_load_matrix_to_LP(Prob_trip, LP));
    // convert matrix to a CSC matrix
    OK(SPEX_matrix_copy(&Prob, SPEX_CSC, SPEX_MPZ, Prob_trip, option));
    // free the triplet matrix
    OK(SPEX_matrix_free(&Prob_trip, option));
    n = Prob->m;

    // read arrays b, lb, ub and c to LP
    SPEX_finish_LP(LP, file_name, option);

    //--------------------------------------------------------------------------
    // generate initial basis matrix
    //--------------------------------------------------------------------------
    glp_adv_basis(LP, 0);
    glp_factorize(LP); // in order to use glp_get_bhead
    basis = (int64_t*) SPEX_malloc(n * sizeof(int64_t));
    if (!basis)
    {
        FREE_WORKSPACE;
        return 0;
    }

    int sgn;
    mpz_t tmpz;
    OK(SPEX_mpz_init(tmpz));
    OK(SPEX_mpq_get_den(tmpz, Prob->scale));
    OK(SPEX_mpz_cmp_ui(&sgn, tmpz, 1));
    mpz_clear(tmpz);
    if (sgn != 0)
    {
        OK(SPEX_gmp_printf("scale is %Qd, which is not 1\n", Prob->scale));
        FREE_WORKSPACE;
        return 0;
    }

    // allocate A2 with n sparse vectors with initially 0 nnz
    OK(SPEX_mat_alloc(&A2, n, n, true));
    for (i = 0; i < n; i++)
    {
        basis[i] = glp_get_bhead(LP, i+1);
        if (basis[i] <= n)
        {
            OK(SPEX_vector_realloc(A2->v[i], 1));
            mpq_get_num(A2->v[i]->x[0], Prob->scale);
            A2->v[i]->i[0] = basis[i]-1;
            A2->v[i]->nz = 1;
        }
        else
        {
            j = basis[i]-n-1;
            OK(SPEX_vector_realloc(A2->v[i], Prob->p[j+1]-Prob->p[j]));
            nz = 0;
            for (p = Prob->p[j]; p < Prob->p[j+1]; p++)
            {
                A2->v[i]->i[nz] = Prob->i[p];
                OK(SPEX_mpz_set(A2->v[i]->x[nz], SPEX_1D(Prob, p, mpz)));
                nz++;
            }
            A2->v[i]->nz = nz;
        }
    }
    OK(SPEX_mat_to_CSC(&A1, A2, NULL, true, option));

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
    S  = SPEX_create_mpq_array(3*n);
    if (!d || !sd || !S || !P || !Q || !P_inv || !Q_inv)
    {
        FREE_WORKSPACE;
        return 0;
    }
    for (i = 0; i < n; i++)
    {
        Q[i] = analysis->q[i];
        P_inv[i] = P1_inv[i];
        P[P_inv[i]] = i;
        Q_inv[Q[i]] = i;
        OK(SPEX_mpz_set(sd[i], SPEX_1D(rhos, i, mpz)));
        OK(SPEX_mpz_set(d[i], sd[i]));
        OK(SPEX_mpq_set_ui(SPEX_LUU_2D(S, 1, i), 1, 1));
        OK(SPEX_mpq_set_ui(SPEX_LUU_2D(S, 2, i), 1, 1));
        OK(SPEX_mpq_set_ui(SPEX_LUU_2D(S, 3, i), 1, 1));
    }

    // convert the factorization to SPEX_mat to be used in the update process
    OK(SPEX_CSC_to_mat(&L2, P,    true,  L1, option));
    OK(SPEX_CSC_to_mat(&U2, Q,    false, U1, option));// U2 stored in row-wise

    // allocate space for vk
    OK(SPEX_vector_alloc(&vk, 0, true));

    int64_t k = 0, new_col = 0;
    while (glp_get_status(LP) != GLP_OPT)
    {
    GOTCHA;
        //----------------------------------------------------------------------
        // run one iteration of simplex and find the new basis
        //----------------------------------------------------------------------
        glp_simplex(LP, &parm);
        k = -1;
        for (i = 0; i < n; i++)
        {
            new_col = glp_get_bhead(LP, i+1);
            if (basis[i] != new_col)
            {
                if (k > -1)
                {
                    printf("replacing more than 1 col\n");
                    FREE_WORKSPACE;
                    return 0;
                }
                k = i;
                basis[i] = new_col;
            }
        }
        if (k == -1)
        {
            printf("same basis matrix\n");
            //continue;
            FREE_WORKSPACE;
            return 0;
        }

        //----------------------------------------------------------------------
        // construct vk
        //----------------------------------------------------------------------
        new_col = basis[k];
        if (new_col <= n)
        {
            vk->i[0] = new_col-1;
            mpq_get_num(vk->x[0], Prob->scale);
            vk->nz = 1;
        }
        else
        {
            new_col = new_col-n-1;
            vk->nz = Prob->p[new_col+1]-Prob->p[new_col];
            if (vk->nzmax < vk->nz)
            {
                OK(SPEX_vector_realloc(vk, vk->nz));
            }
            i = 0;
            for (p = Prob->p[new_col]; p < Prob->p[new_col+1]; p++)
            {
                vk->i[i] = Prob->i[p];
                OK(SPEX_mpz_set(vk->x[i], SPEX_1D(Prob, p, mpz)));
                i++;
            }
            vk->nz = i;
        }
        /*int num[5]={0,0,0,0,0};
        for (int i = 1; i <= n; i++)
        {
            int stat=glp_get_row_stat(LP,i);
            num[stat-1]++;
        }
        for (int i = 1; i <= Prob->n; i++)
        {
            int stat=glp_get_col_stat(LP,i);
            num[stat-1]++;
        }
        printf("bs:%d nl:%d nu:%d nf:%d ns:%d\n",num[0],num[1],num[2],num[3],num[4]);*/

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

        end_llu = clock();

        //----------------------------------------------------------------------
        // perform LU update for matrix A2->A1
        //----------------------------------------------------------------------
        start_luu = clock();

        OK(SPEX_LUU(A2, L2, U2, d, sd, S, P, P_inv, Q, Q_inv, &vk, k, option));

        end_luu = clock();

        //----------------------------------------------------------------------
        // print results
        //----------------------------------------------------------------------
        // Timing stats
        double t_llu = (double) (end_llu-start_llu)/CLOCKS_PER_SEC;
        double t_luu = (double) (end_luu - start_luu) / CLOCKS_PER_SEC;
        printf("\nSPEX Left LU Factorization time: \t%lf", t_llu);
        printf("\nSPEX LU Update time: \t\t\t%lf\n\n", t_luu);
    }

    FREE_WORKSPACE;
    return 0;
}

