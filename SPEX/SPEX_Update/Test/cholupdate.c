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
    /*fclose(out_file);*/                            \
    if (mat_file!=NULL) fclose(mat_file);        \
    SPEX_FREE(option);                           \
    SPEX_finalize() ;                            \
}

#include "test.h"
#include "SPEX_Chol.h"
#include <assert.h>

int main( int argc, char* argv[])
{
    //char *prob_name = "lp_80bau3b";
    char *prob_name = "lp_25fv47";
    //char *prob_name = "lp_afiro"; // optimal: -4.6475314E+02
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
    int64_t n, i, j, p, n_col, target = -1, target_sym = -1, i1, j1, pi, pj;
    int sgn;
    SPEX_options *option = NULL;
    SPEX_Chol_analysis *analysis = NULL;
    SPEX_matrix *Prob_A = NULL, *Prob_b = NULL, *Prob_c = NULL;
    SPEX_matrix *L1 = NULL, *rhos1 = NULL, *A1 = NULL, *tmpA = NULL;
    SPEX_matrix *rhos2 = NULL;
    SPEX_matrix *L2 = NULL, *A2 = NULL;
    SPEX_vector *w = NULL;
    mpz_t tmpz, tmpz1;
    int64_t *P2_inv = NULL, *P2 = NULL;
    int64_t *basis = NULL, *used= NULL;
    clock_t start1, start2, end1, end2;
    FILE *mat_file = NULL;
    glp_prob *LP;
    glp_smcp parm;
    LP = glp_create_prob();
    glp_init_smcp(&parm);
    parm.it_lim = 1;
    OK(SPEX_create_default_options(&option));
    OK(SPEX_mpz_init(tmpz));
    OK(SPEX_mpz_init(tmpz1));

    //--------------------------------------------------------------------------
    // read in matrix
    //--------------------------------------------------------------------------
    char file_name[1000] = "TestMats/";
    strcat(file_name, prob_name);
    // read matrix A, arrays b, lb, ub and c to LP
    SPEX_construct_LP(LP, &Prob_A, &Prob_b, &Prob_c, file_name, option);
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

    OK(SPEX_mpq_get_den(tmpz, Prob_A->scale));
    OK(SPEX_mpz_cmp_ui(&sgn, tmpz, 1));
    if (sgn != 0)
    {
        OK(SPEX_gmp_printf("scale is %Qd, whose den is not 1\n",
            Prob_A->scale));
        FREE_WORKSPACE;
        return 0;
    }

    // allocate A2 with n sparse vectors with initially n nnz
    OK(SPEX_matrix_allocate(&A2, SPEX_DYNAMIC_CSC, SPEX_MPZ, n, n, 0, false,
        true, option));
    for (i = 0; i < n; i++)
    {
        OK(SPEX_vector_realloc(A2->v[i], n, option));
    }
    printf("set of basic variables found, now computing A=B*B^T....\n");
    for (i = 0; i < n; i++)
    {
        if (basis[i] < n)
        {
            j = basis[i];
            target = -1;
            for (p = 0; p < A2->v[j]->nz; p++)
            {
                if (A2->v[j]->i[p] == j) {target = p; break;}
            }
            if (target == -1)
            {
                target = A2->v[j]->nz;
                A2->v[j]->i[target] = j;
                A2->v[j]->nz++;
            }
            mpq_get_num(tmpz, Prob_A->scale);
            OK(SPEX_mpz_addmul(A2->v[j]->x[target], tmpz, tmpz));
        }
        else
        {
            j = basis[i]-n;
            used[j] = 1;

            for (pj = Prob_A->p[j]; pj < Prob_A->p[j+1]; pj++)
            {
                j1 = Prob_A->i[pj];

                target = -1;
                for (p = 0; p < A2->v[j1]->nz; p++)
                {
                    if (A2->v[j1]->i[p] == j1) {target = p; break;}
                }
                if (target == -1)
                {
                    target = A2->v[j1]->nz;
                    A2->v[j1]->i[target] = j1;
                    A2->v[j1]->nz++;
                }
                OK(SPEX_mpz_addmul(A2->v[j1]->x[target],
                                   SPEX_1D(Prob_A, pj, mpz),
                                   SPEX_1D(Prob_A, pj, mpz)));

                for (pi = pj+1; pi < Prob_A->p[j+1]; pi++)
                {
                    i1 = Prob_A->i[pi];

                    target = -1;
                    for (p = 0; p < A2->v[j1]->nz; p++)
                    {
                        if (A2->v[j1]->i[p] == i1) {target = p; break;}
                    }
                    if (target == -1)
                    {
                        target = A2->v[j1]->nz;
                        A2->v[j1]->i[target] = i1;
                        A2->v[j1]->nz++;
                        target_sym = A2->v[i1]->nz;
                        A2->v[i1]->i[target_sym] = j1;
                        A2->v[i1]->nz++;
                    }
                    else
                    {
                        target_sym = -1;
                        for (p = 0; p < A2->v[i1]->nz; p++)
                        {
                            if (A2->v[i1]->i[p] == j1) {target_sym = p; break;}
                        }
                        if (target_sym == -1)
                        {
                            printf("this shouldn't happen, ");
                            printf("have fun finding the bug\n");
                            return 0;
                        }
                    }
                    OK(SPEX_mpz_addmul(A2->v[j1]->x[target],
                                       SPEX_1D(Prob_A, pj, mpz),
                                       SPEX_1D(Prob_A, pi, mpz)));
                    OK(SPEX_mpz_set(A2->v[i1]->x[target_sym],
                                    A2->v[j1]->x[target]));
                }
            }
        }
    }
    OK(SPEX_mpq_mul(A2->scale, Prob_A->scale, Prob_A->scale));
    OK(SPEX_matrix_copy(&A1, SPEX_CSC, SPEX_MPZ, A2, option));

    //--------------------------------------------------------------------------
    // perform Cholesky factorization
    //--------------------------------------------------------------------------
    printf("compute initial Cholesky factorization for A....\n");
    OK(SPEX_Chol_preorder(&analysis, A1, option));
    OK(SPEX_Chol_permute_A(&tmpA, A1, analysis));
    bool left_looking = true;// True = left, false = up
    OK(SPEX_Chol_Factor(&L1, &rhos1, analysis, tmpA, left_looking, option));

    //--------------------------------------------------------------------------
    // generate initial inputs for LU update
    //--------------------------------------------------------------------------
    printf("generating inputs for Cholesky rank-1 update....\n");
    // generate permutation vectors P, Q, P_inv, Q_inv and vectors sd
    P2     = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    P2_inv = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    if (!P2 || !P2_inv)
    {
        FREE_WORKSPACE;
        return 0;
    }
    OK(SPEX_matrix_allocate(&rhos2, SPEX_DENSE, SPEX_MPZ, n, 1, n, false, true,
        option));
    for (i = 0; i < n; i++)
    {
        P2[i] = analysis->q[i];
        P2_inv[P2[i]] = i;
        OK(SPEX_mpz_set(SPEX_1D(rhos2, i, mpz), SPEX_1D(rhos1, i, mpz)));
    }

    // convert the factorization to SPEX_mat to be used in the update process
    OK(SPEX_matrix_copy(&L2, SPEX_DYNAMIC_CSC, SPEX_MPZ, L1, option));
    OK(SPEX_permute_row(L2, P2, option));
    OK(SPEX_matrix_canonicalize(L2, P2, option));

    // allocate space for scattered vector w
    OK(SPEX_vector_allocate(&w, 0, true, option));
    OK(SPEX_vector_realloc(w, n, option));

    //--------------------------------------------------------------------------
    // perform update
    //--------------------------------------------------------------------------
    int64_t new_col = -1;
    for (i = 0; i < n_col; i++)
    {
        if (used[i] != 1)
        {
            new_col = i;
            j = 0;
            for (p = Prob_A->p[new_col]; p < Prob_A->p[new_col+1]; p++)
            {
#if 0
                for (; j < Prob_A->i[p]; j++)
                {
                    OK(SPEX_mpz_set_ui(y->v[0]->x[j], 0));
                }
                OK(SPEX_mpz_set(w->x[j], SPEX_1D(Prob_A, p, mpz)));// dense w
                j++;
#else
                j = Prob_A->i[p];
                OK(SPEX_mpz_set(w->x[j], SPEX_1D(Prob_A, p, mpz)));// dense w
#endif
            }
            used[i] = 1;
            break;
        }

    }
    printf("computing Cholesky rank-1 update if column %ld is add to B...\n",
        new_col);
    start1 = clock();
    OK(SPEX_Update_Chol_Rank1(L2, rhos2, P2, P2_inv, w, 1));
    end1 = clock();

    //--------------------------------------------------------------------------
    // compute updated Cholesky using direct factorization
    //--------------------------------------------------------------------------
    printf("computing updated Cholesky using direct factorization...\n");
    j = new_col;
    for (pj = Prob_A->p[j]; pj < Prob_A->p[j+1]; pj++)
    {
        j1 = Prob_A->i[pj];

        target = -1;
        for (p = 0; p < A2->v[j1]->nz; p++)
        {
            if (A2->v[j1]->i[p] == j1) {target = p; break;}
        }
        if (target == -1)
        {
            target = A2->v[j1]->nz;
            A2->v[j1]->i[target] = j1;
            A2->v[j1]->nz++;
        }
        OK(SPEX_mpz_addmul(A2->v[j1]->x[target],
                           SPEX_1D(Prob_A, pj, mpz),
                           SPEX_1D(Prob_A, pj, mpz)));

        for (pi = pj+1; pi < Prob_A->p[j+1]; pi++)
        {
            i1 = Prob_A->i[pi];

            target = -1;
            for (p = 0; p < A2->v[j1]->nz; p++)
            {
                if (A2->v[j1]->i[p] == i1) {target = p; break;}
            }
            if (target == -1)
            {
                target = A2->v[j1]->nz;
                A2->v[j1]->i[target] = i1;
                A2->v[j1]->nz++;
                target_sym = A2->v[i1]->nz;
                A2->v[i1]->i[target_sym] = j1;
                A2->v[i1]->nz++;
            }
            else
            {
                target_sym = -1;
                for (p = 0; p < A2->v[i1]->nz; p++)
                {
                    if (A2->v[i1]->i[p] == j1) {target_sym = p; break;}
                }
                if (target_sym == -1)
                {
                    printf("this shouldn't happen, ");
                    printf("have fun finding the bug\n");
                    return 0;
                }
            }
            OK(SPEX_mpz_addmul(A2->v[j1]->x[target],
                               SPEX_1D(Prob_A, pj, mpz),
                               SPEX_1D(Prob_A, pi, mpz)));
            OK(SPEX_mpz_set(A2->v[i1]->x[target_sym],
                            A2->v[j1]->x[target]));
        }
    }
    SPEX_matrix_free(&A1, NULL);
    OK(SPEX_matrix_copy(&A1, SPEX_CSC, SPEX_MPZ, A2, option));
    SPEX_FREE(analysis);
    SPEX_matrix_free(&tmpA, NULL);
    SPEX_matrix_free(&rhos1, NULL);

    // perform Cholesky factorization
    start2 = clock();
    OK(SPEX_Chol_preorder(&analysis, A1, option));
    OK(SPEX_Chol_permute_A(&tmpA, A1, analysis));
    OK(SPEX_Chol_Factor(&L1, &rhos1, analysis, tmpA, left_looking, option));
    end2 = clock();

    for (j = 0; j < n; j++)
    {
        for (p = L1->p[j]; p < L1->p[j+1]; p++)
        {
            i = L1->i[p];
            //gmp_printf("%ld %Zd\n",i, L1->x.mpz[p]);
            OK(SPEX_mpz_sgn(&sgn, L1->x.mpz[p]));
            if (sgn == 0) {continue;}

            bool found = false;
            for (pi = 0; pi < L2->v[j]->nz; pi++)
            {
                //printf("P[%ld]=%ld %ld\n",i, P2[i],L2->v[j]->i[pi]);
                if (L2->v[j]->i[pi] == P2[i])
                {
                    mpq_get_num(tmpz1, L2->v[j]->scale);
                    OK(SPEX_mpz_mul(tmpz, L2->v[j]->x[pi], tmpz1));
                    OK(SPEX_mpq_get_den(tmpz1, L2->v[j]->scale));
                    OK(SPEX_mpz_divexact(tmpz, tmpz, tmpz1));
                    OK(SPEX_mpz_cmp(&sgn, tmpz, L1->x.mpz[p]));
                    if (sgn != 0)
                    {
                        printf("found different value for L(%ld, %ld)\n",
                            i, j);
                        return 0;
                    }
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                printf("failed to find corresponding L(%ld, %ld)\n", i, j);
                return 0;
            }
        }
    }
    //----------------------------------------------------------------------
    // print results
    //----------------------------------------------------------------------
    // Timing stats
    double t1 = (double) (end1 - start1)/CLOCKS_PER_SEC;
    double t2 = (double) (end2 - start2) / CLOCKS_PER_SEC;

    printf("test succeeded!!\n");
    printf("\nSPEX Cholesky Factorization time: \t%lf", t2);
    printf("\nSPEX Cholesky rank-1 Update time: \t%lf\n", t1);


    FREE_WORKSPACE;
    return 0;
}
