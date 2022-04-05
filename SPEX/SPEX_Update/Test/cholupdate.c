//------------------------------------------------------------------------------
// SPEX_Update/Test/cholupdate.c: test Cholesky rank-1 update functionality
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
    glp_delete_prob(LP);                         \
    glp_free_env();                              \
    if (result_file != NULL) fclose(result_file);\
    SPEX_matrix_free(&Prob_A, option);           \
    SPEX_matrix_free(&Prob_b, option);           \
    SPEX_matrix_free(&Prob_c, option);           \
    SPEX_matrix_free(&A1, option);               \
    SPEX_matrix_free(&A2, option);               \
    SPEX_matrix_free(&A3, option);               \
    SPEX_matrix_free(&A0, option);               \
    SPEX_matrix_free(&w, option);                \
    SPEX_symbolic_analysis_free(&analysis, option);\
    SPEX_factorization_free(&F1, option);        \
    SPEX_factorization_free(&F2, option);        \
    SPEX_factorization_free(&F3, option);        \
    SPEX_factorization_free(&F_update, option);  \
    SPEX_factorization_free(&Ftmp, option);      \
    SPEX_FREE(option);                           \
    mpz_clear(tmpz);                             \
    SPEX_FREE(basis);                            \
    SPEX_FREE(used);                             \
    SPEX_finalize() ;                            \
}
#define PRINT_TO_FILE

#include "test.h"
#include "SPEX_Chol.h"
#include <assert.h>

int main( int argc, char* argv[])
{
    //char *prob_name = "lp_80bau3b";
    //char *prob_name = "lp_25fv47";
    char *prob_name = "lp_afiro"; // optimal: -4.6475314E+02
    //char *prob_name = "aa5";
    if (argc >= 2)
    {
        prob_name = argv[1];
    }
    SPEX_info info;
    //------------------------------------------------------------------
    // Initialize for SPEX library
    //------------------------------------------------------------------

    SPEX_initialize () ;

    //------------------------------------------------------------------
    // Initialize and allocate workspace
    //------------------------------------------------------------------
    int64_t n, i, j, p, n_col, target = -1, sigma;
    int sgn;
    SPEX_options *option = NULL;
    SPEX_symbolic_analysis *analysis = NULL;
    SPEX_factorization *F1 = NULL, *F2 = NULL, *F3 = NULL, *F_update = NULL,
                *Ftmp = NULL;
    SPEX_matrix *Prob_A = NULL, *Prob_b = NULL, *Prob_c = NULL;
    SPEX_matrix *A1 = NULL, *A2 = NULL, *A3 = NULL;
    SPEX_matrix *A0 = NULL, *w = NULL;
    mpz_t tmpz;
    int64_t *basis = NULL, *used= NULL;
    clock_t start, end, start1, start2, end1, end2;
    char file_name[1000];
    FILE *result_file = NULL;
    glp_prob *LP;
    glp_smcp parm;
    LP = glp_create_prob();
    glp_init_smcp(&parm);
    parm.it_lim = 1;
    OK(SPEX_create_default_options(&option));
    OK(SPEX_mpz_init(tmpz));

    //--------------------------------------------------------------------------
    // open output file
    //--------------------------------------------------------------------------
#ifdef PRINT_TO_FILE
    sprintf(file_name, "Results/LPnetlib_CholUpdate/%s_chol.txt",prob_name);
    result_file = fopen(file_name, "w");
   if (result_file == NULL)
    {
        printf("Could not open file: %s!\n", file_name);
        FREE_WORKSPACE;
        return 0;
    }
#endif

    //--------------------------------------------------------------------------
    // read in matrix
    //--------------------------------------------------------------------------
    sprintf(file_name, "TestMats/LPnetlib/%s/%s",prob_name,prob_name);
    // read matrix A, arrays b, lb, ub and c to LP
    double z0 = 0;
    SPEX_construct_LP(LP, &Prob_A, &Prob_b, &Prob_c, &z0, file_name, option);
    n = Prob_A->m;
    n_col = Prob_A->n; // number of columns in matrix A

    //--------------------------------------------------------------------------
    // generate initial basis matrix
    //--------------------------------------------------------------------------
    basis = (int64_t*) SPEX_malloc(n_col * sizeof(int64_t));
    used  = (int64_t*) SPEX_calloc(n_col,  sizeof(int64_t));
    if (!basis || !used)
    {
        FREE_WORKSPACE;
        return 0;
    }

    printf("getting initial basic variables....\n");
    glp_adv_basis(LP, 0);
    glp_factorize(LP); // in order to use glp_get_bhead
    parm.msg_lev = GLP_MSG_OFF;
    while (glp_get_status(LP) != GLP_FEAS)
    {
        glp_simplex(LP, &parm);
    }
    if (glp_get_status(LP) == GLP_OPT)
    {
        printf("LP is optimized!\n");
        FREE_WORKSPACE;
        return 0;
    }
    for (i = 0; i < n; i++)
    {
        basis[i] = glp_get_bhead(LP, i+1)-1;
#if 0
        if (basis[i] < n && glp_get_mat_row(LP, basis[i]+1, NULL, NULL) != 0)
        {
            //printf("basis[%ld]=%ld\n",i,basis[i]);
            GOTCHA;
            return 0;;
        }
#endif
    }

    OK(SPEX_mpq_get_den(tmpz, Prob_A->scale));
    OK(SPEX_mpz_cmp_ui(&sgn, tmpz, 1));
    if (sgn != 0)
    {
        OK(SPEX_gmp_printf("scale is %Qd, whose den is not 1\n",
            Prob_A->scale));
        FREE_WORKSPACE;
        return 0;
    }

    // allocate A0 with n sparse vectors with initially n nnz
    OK(SPEX_matrix_allocate(&A0, SPEX_DYNAMIC_CSC, SPEX_MPZ, n, n, 0, false,
        true, option));
    for (i = 0; i < n; i++)
    {
        OK(SPEX_vector_realloc(A0->v[i], n, option));
    }

    //--------------------------------------------------------------------------
    // compute A1 as B*B^T 
    //--------------------------------------------------------------------------
    printf("set of basic variables found, now computing A=B*B^T....\n");
    start = clock();
    for (i = 0; i < n; i++)
    {
        if (basis[i] < n)
        {
            j = basis[i];
            target = -1;
            for (p = 0; p < A0->v[j]->nz; p++)
            {
                if (A0->v[j]->i[p] == j) {target = p; break;}
            }
            if (target == -1)
            {
                target = A0->v[j]->nz;
                A0->v[j]->i[target] = j;
                A0->v[j]->nz++;
            }
            mpq_get_num(tmpz, Prob_A->scale);
            OK(SPEX_mpz_addmul(A0->v[j]->x[target], tmpz, tmpz));
        }
        else
        {
            j = basis[i]-n;
            used[j] = 1;

            OK(SPEX_A_plus_vvT(A0, Prob_A, j));
        }
    }
    OK(SPEX_mpq_mul(A0->scale, Prob_A->scale, Prob_A->scale));
    OK(SPEX_matrix_copy(&A1, SPEX_CSC, SPEX_MPZ, A0, option));
    GOTCHA;

    //--------------------------------------------------------------------------
    // find the most sparse and dense vectors from remaining cols
    //--------------------------------------------------------------------------
    int64_t sparsest = -1, sparse_nz = n, densest = -1, dense_nz = 0;
    for (i = 0; i < n_col; i++)
    {
        if (used[i] != 1)
        {
            int64_t Ai_nz = Prob_A->p[i+1] - Prob_A->p[i];
            if (Ai_nz < sparse_nz)
            {
                sparsest = i;
                sparse_nz = Ai_nz;
            }
            if (Ai_nz > dense_nz)
            {
                densest = i;
                dense_nz = Ai_nz;
            }
        }

    }

    if (sparsest == -1 || densest == -1)
    {
        printf("cannot find any unused column\n");
        return 0;
    }

    //--------------------------------------------------------------------------
    // compute A2 as A1+v_sparsest*v_sparsest^T
    //--------------------------------------------------------------------------
    OK(SPEX_A_plus_vvT(A0, Prob_A, sparsest));
    OK(SPEX_matrix_copy(&A2, SPEX_CSC, SPEX_MPZ, A0, option));

    //--------------------------------------------------------------------------
    // compute A3 as A2+v_densest*v_densest^T
    //--------------------------------------------------------------------------
    OK(SPEX_A_plus_vvT(A0, Prob_A, densest));
    OK(SPEX_matrix_copy(&A3, SPEX_CSC, SPEX_MPZ, A0, option));
    end = clock();
    printf("time to compute BB^T: %lf\n", (double) (end - start)/CLOCKS_PER_SEC);

    //--------------------------------------------------------------------------
    // perform Cholesky factorization for A1, A2, A3 with same pre-ordering
    //--------------------------------------------------------------------------
    printf("compute initial Cholesky factorization for A....\n");
    option->algo = SPEX_CHOL_LEFT;// or SPEX_CHOL_UP
    GOTCHA;
    OK(SPEX_Chol_analyze(&analysis, A1, option));
    GOTCHA;
    OK(SPEX_Chol_factorize(&F1, analysis, A1, option));
    GOTCHA;
    OK(SPEX_Chol_factorize(&F2, analysis, A2, option));
    GOTCHA;
    OK(SPEX_Chol_factorize(&F3, analysis, A3, option));
    GOTCHA;

    //--------------------------------------------------------------------------
    // generate initial inputs for LU update
    //--------------------------------------------------------------------------
    printf("generating inputs for Cholesky rank-1 update....\n");
    // make a copy of F1 as F_update without converting to updatable
    F_update = (SPEX_factorization*) SPEX_calloc(1, sizeof(SPEX_factorization));
    if (F_update == NULL)
    {
        FREE_WORKSPACE  ;
        return 0;
    }
    // set factorization kind
    F_update->kind = SPEX_CHOLESKY_FACTORIZATION;
    // Allocate and set scale_for_A
    OK(SPEX_mpq_init(F_update->scale_for_A));
    OK(SPEX_mpq_set (F_update->scale_for_A, F1->scale_for_A));

    // allocate Pinv and P
    F_update->Pinv_perm = (int64_t*) SPEX_malloc (n * sizeof(int64_t));
    F_update->P_perm =    (int64_t*) SPEX_malloc (n * sizeof(int64_t));

    if (!(F_update->Pinv_perm) || !(F_update->P_perm))
    {
        // out of memory: free everything and return
        FREE_WORKSPACE  ;
        return 0;
    }

    // copy column permutation, rhos and L from F1
    memcpy(F_update->P_perm,    F1->P_perm,    n * sizeof(int64_t));
    memcpy(F_update->Pinv_perm, F1->Pinv_perm, n * sizeof(int64_t));
    OK(SPEX_matrix_copy(&F_update->rhos, SPEX_DENSE, SPEX_MPZ, F1->rhos,
        option));
    OK(SPEX_matrix_copy(&F_update->L, SPEX_CSC, SPEX_MPZ, F1->L, option));

    // allocate space for scattered vector w
    OK(SPEX_matrix_allocate(&w, SPEX_DYNAMIC_CSC, SPEX_MPZ, n, 1, 0, false,
        true, option));
    OK(SPEX_vector_realloc(w->v[0], n, option));

    for (int iter = 0; iter < 4; iter++)
    {
    GOTCHA;
        if (iter < 2) sigma = 1; // the first two iterations are updates
        else sigma = -1;// while the last two iterations are downdates

        //----------------------------------------------------------------------
        // get the vector w
        //----------------------------------------------------------------------
        int64_t new_col;
        if (iter == 1 || iter == 2) new_col = densest;
        else new_col = sparsest;

        j = 0;
        for (p = Prob_A->p[new_col]; p < Prob_A->p[new_col+1]; p++)
        {
            OK(SPEX_mpz_set(w->v[0]->x[j], SPEX_1D(Prob_A, p, mpz)));
            w->v[0]->i[j] = Prob_A->i[p];
            j++;
        }
        w->v[0]->nz = j;
        // used[new_col] = 1;

        //----------------------------------------------------------------------
        // perform update
        //----------------------------------------------------------------------
        printf("computing Cholesky rank-1 update if col %ld is %s B...\n",
            new_col, iter < 2?"added to":"removed from");
        start1 = clock();
        OK(SPEX_Update_Chol_Rank1(F_update, w, sigma, option));
        end1 = clock();

        //----------------------------------------------------------------------
        // compute updated Cholesky using direct factorization
        //----------------------------------------------------------------------
        printf("computing updated Cholesky using direct factorization...\n");
        SPEX_symbolic_analysis_free(&analysis, option);
        SPEX_factorization_free(&Ftmp, option);

        // perform Cholesky factorization, skip iter == 2 which is same as 0
        start2 = clock();
        if (iter == 0)
        {
            OK(SPEX_Chol_analyze(&analysis, A2, option));
            OK(SPEX_Chol_factorize(&Ftmp, analysis, A2, option));
        }
        else if (iter == 1)
        {
            OK(SPEX_Chol_analyze(&analysis, A3, option));
            OK(SPEX_Chol_factorize(&Ftmp, analysis, A3, option));
        }
        else if(iter == 3)
        {
            OK(SPEX_Chol_analyze(&analysis, A1, option));
            OK(SPEX_Chol_factorize(&Ftmp, analysis, A1, option));
        }
        end2 = clock();

        bool Isequal;
        size_t s, L_sum_size = 0, L_update_sum_size = 0;
        int64_t L_nnz, L_update_nnz = 0;
        SPEX_matrix *L;
        if (iter == 0 || iter == 2)
        {
            L = F2->L;
        }
        else if (iter == 1)
        {
            L = F3->L;
        }
        else // if(iter == 3)
        {
            L = F1->L;
        }
        SPEX_matrix_equal(&Isequal, L, F_update->L, F_update->P_perm);
        if (!Isequal) return 0;

        for (p = 0; p < L->p[n]; p++)
        {
            OK(SPEX_mpz_sizeinbase(&s, L->x.mpz[p], 2));
            L_sum_size += s;
        }
        L_nnz = L->p[n];

        for (i = 0; i < n; i++)
        {
            L_update_nnz += F_update->L->v[i]->nz;
            for (p = 0; p < F_update->L->v[i]->nz; p++)
            {
                OK(SPEX_mpz_sizeinbase(&s, F_update->L->v[i]->x[p], 2));
                L_update_sum_size += s;
            }
        }
        //----------------------------------------------------------------------
        // print results
        //----------------------------------------------------------------------
        // Timing stats
        double t1 = (double) (end1 - start1)/CLOCKS_PER_SEC;
        double t2 = (double) (end2 - start2) / CLOCKS_PER_SEC;

        printf("test succeeded!!");
        printf("\n\t\t time \t\tL_nnz \tsum(Lb)");
        printf("\nFactorization: \t%lf \t%ld \t%ld", t2,L_nnz,L_sum_size);
        printf("\nrank-1 Update: \t%lf \t%ld \t%ld\n\n", t1,L_update_nnz,
            L_update_sum_size);
#ifdef PRINT_TO_FILE
        if (iter == 2)
        {
            fprintf(result_file,"\t%lf \t%ld \t%ld",t1,L_update_nnz,
                L_update_sum_size);
        }
        else
        {
            fprintf(result_file,"\t%lf \t%ld \t%ld",t2,L_nnz,L_sum_size);
            fprintf(result_file,"\t%lf \t%ld \t%ld",t1,L_update_nnz,
                L_update_sum_size);
        }
#endif
    }

    FREE_WORKSPACE;
    return 0;
}
