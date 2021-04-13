//------------------------------------------------------------------------------
// SPEX_CHOLMOD/Tcov/tcov_test.c: test coverage for SPEX_CHOLMOD
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/*
 * test coverage
 * The format of the input file should be as follow:
 * First line of the file give the size of the frame matrix (n), the number of
 * its nnz, the index of the column (k) that is going to be replaced (k
 * will be 0-based), and the number of nnz in the new column.
 * Then, the file should list all the entries of the frame matrix in triplet
 * form .
 * Finally, the entries of the new column are listed with its row index and
 * corresponding value
 * The frame matrix contains both L and U. This code requires all of its
 * entries to be either 0 or 1 and its diagonal must be all 1s.
 */

#define SPEX_FREE_ALL                            \
{                                                \
    if (mat_file != NULL) fclose(mat_file);      \
    SPEX_FREE(option);                           \
    SPEX_matrix_free(&L);                        \
    SPEX_matrix_free(&U);                        \
    SPEX_matrix_free(&A);                        \
    SPEX_vector_free(&vk);                       \
    spex_delete_mpz_array(&d, An);                \
    spex_delete_mpz_array(&sd, An);               \
    spex_delete_mpq_array(&S, 3*An);              \
    SPEX_FREE(P);                                \
    SPEX_FREE(P_inv);                            \
    SPEX_FREE(Q);                                \
    SPEX_FREE(Q_inv);                            \
    SPEX_MPZ_CLEAR(tmpz);                        \
    SPEX_MPQ_CLEAR(tmpq);                        \
    SPEX_finalize() ;                            \
}


#include "tcov_malloc_test.h"

#define TEST_CHECK(method)                       \
{                                                \
    info = (method) ;                            \
    if (info != SPEX_OK)                         \
    {                                            \
        SPEX_PRINT_INFO (info) ;                 \
        SPEX_FREE_ALL;printf("n=%ld\n",n);                           \
        return 0 ;/*continue; */                               \
    }                                            \
}

#define TEST_CHECK_FAILURE(method)               \
{                                                \
    info = (method) ;                            \
    if (info != SPEX_INCORRECT_INPUT && info != SPEX_SINGULAR) \
    {                                            \
        SPEX_PRINT_INFO (info) ;                 \
        SPEX_FREE_ALL ;                          \
        continue ;                               \
    }                                            \
    else                                         \
    {                                            \
        printf("Expected failure at line %d\n", __LINE__);\
    }                                            \
}

#define MAX_MALLOC_COUNT 1000

int64_t Aold[16] = {3, 5, 6, 7, 11, 2, -7, 10, 8, 3, -2, -2, 7, 5, 1, -6};//CSC
int64_t Lold[10] = {3, 5, 6, 7, -49, -87, -47, -17, 527, 884};// CSC
int64_t Uold[10] = {3, 11, 8, 7, -49, -31, -20, -17, 57,884};//CSR
int64_t Ak_new[4]= {1,4,7,11}; // new column
#include <assert.h>

int main( int argc, char* argv[])
{
    SPEX_info info;
    int sgn;
    int64_t n, n_min = 4, n_max, An=0, Anz, vk_nz, i, j, k = 0, p;
    bool read_matrix = true;
    GOTCHA;

    if (argc >= 2)
    {
        n_min = atoi(argv[1]);
    }
    n_max = n_min+1;
    if (argc >= 3)
    {
        n_max = atoi(argv[2]);
    }
    if (n_min == 0)
    {
        read_matrix = true;
        if (argc > 3 && argc < n_max+3) {return 0;}
    }
    else
    {
        read_matrix = false;
    }

    for (n = n_min; n < n_max; n++)
    {
        GOTCHA;
        for (int64_t loop = 0; loop < (read_matrix?1:100); loop++)
        {
            printf("\n\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n+++++++++++++++++++++++++++++++++case %ld+++++++++++++++++++++++++++++++++\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n\n",loop);
            //------------------------------------------------------------------
            // Initialize SPEX CHOLMOD process
            //------------------------------------------------------------------

            SPEX_initialize_expert (tcov_malloc, tcov_calloc,
                tcov_realloc, tcov_free) ;

            info = SPEX_initialize ( ) ;
            assert (info == SPEX_PANIC) ;
        GOTCHA;

            //------------------------------------------------------------------
            // Allocate memory
            //------------------------------------------------------------------

            SPEX_options* option = SPEX_create_default_options();
            if (!option) return 0;//{continue;}

            SPEX_matrix *L = NULL, *U = NULL, *A = NULL;
            SPEX_vector *vk = NULL;
            mpz_t *d = NULL, *sd = NULL;
            mpq_t *S = NULL;
            int64_t *P = NULL, *P_inv = NULL, *Q = NULL, *Q_inv = NULL;
            mpz_t tmpz; SPEX_MPZ_SET_NULL(tmpz);
            mpq_t tmpq; SPEX_MPQ_SET_NULL(tmpq);
            FILE *mat_file = NULL;
        GOTCHA;

            if (read_matrix)
            {
                // read in L and U
                if (argc > 3)
                {
                    mat_file = fopen(argv[3+n], "r");
                }
                else
                {
                    char mat_name[25];
                    sprintf(mat_name,"Mats4Tcov/mat%ld.txt",n+1);
                    puts(mat_name);
                    mat_file = fopen(mat_name, "r");
                }
                if (mat_file == NULL)
                {
                    perror("Error while opening the file");
                    SPEX_FREE_ALL;
                    return 0;
                }

                // Read in size of matrix & number of nonzeros
                int s = fscanf(mat_file,  "%"PRId64" %"PRId64" %"PRId64""
                    "%"PRId64"\n", &An, &Anz, &k, &vk_nz);
                if (feof(mat_file) || s < 4)
                {
                    printf ("premature end-of-file\n") ;
                    SPEX_FREE_ALL;
                    return 0;
                }

                TEST_CHECK(SPEX_matrix_alloc(&L, An, An, true));
                TEST_CHECK(SPEX_matrix_alloc(&U, An, An, true));
                TEST_CHECK(SPEX_vector_alloc(&vk, An, true));
                for (j = 0; j < An; j++)
                {
                    TEST_CHECK(SPEX_vector_realloc(L->v[j], An-j));
                    TEST_CHECK(SPEX_vector_realloc(U->v[j], An-j));
                }
                // Read in the values from file
                int64_t tmp;
                for (int64_t ai = 0; ai < Anz; ai++)
                {
                    s = fscanf(mat_file, "%"PRId64" %"PRId64" %"PRId64"\n",
                        &i, &j, &tmp);
                    if (feof(mat_file) || s < 3)
                    {
                        printf ("premature end-of-file\n") ;
                        SPEX_FREE_ALL;
                        return 0;
                    }
                    if (j <= i)
                    {
                        TEST_CHECK(SPEX_mpz_set_si(
                            L->v[j]->x[L->v[j]->nz], tmp));
                        L->v[j]->i[L->v[j]->nz] = i;
                        L->v[j]->nz++;
                    }
                    if (j >= i)
                    {
                        TEST_CHECK(SPEX_mpz_set_si(
                            U->v[i]->x[U->v[i]->nz], tmp));
                        U->v[i]->i[U->v[i]->nz] = j;
                        U->v[i]->nz++;
                    }
                }

                for (j = 0; j < vk_nz; j++)
                {
                    s = fscanf(mat_file, "%"PRId64" %"PRId64"\n", &i, &tmp);
                    if ((feof(mat_file) && j != vk_nz-1) || s < 2)
                    {
                        printf ("premature end-of-file\n") ;
                        SPEX_FREE_ALL;
                        return 0;
                    }
                    TEST_CHECK(SPEX_mpz_set_si(vk->x[vk->nz], tmp));
                    vk->i[vk->nz] = i;
                    vk->nz++;
                }
            }
            else
            {
                An = n;
                Anz = 0;
                char *mat_name = "Mats4Tcov/mat.txt";
                mat_file = fopen(mat_name, "w");
                if (mat_file == NULL)
                {
                    perror("Error while opening the file");
                    SPEX_FREE_ALL;
                    return 0;
                }

                // randomly generate L, U and vk
                TEST_CHECK(SPEX_matrix_alloc(&L, An, An, true));
                printf("l0=[\n");
                for (j = 0; j < An; j++)
                {
                    TEST_CHECK(SPEX_vector_realloc(L->v[j], An-j));
                    TEST_CHECK(SPEX_mpz_set_si(L->v[j]->x[0], 1)); // diagnal entry
                    L->v[j]->i[0] = j;
                    p = 1;
                    for (i = 0; i < j; i++){printf("0 ");}
                            printf("1 ");
                    for (i = j+1; i < An; i++)
                    {
                        if (rand() > RAND_MAX/2)
                        {
                            printf("1 ");
                            TEST_CHECK(SPEX_mpz_set_si(L->v[j]->x[p], 1));
                            L->v[j]->i[p] = i;
                            p++;
                        }
                        else
                        {
                            printf("0 ");
                            /*if (rand() > RAND_MAX/10)
                            {
                                TEST_CHECK(SPEX_mpz_set_si(L->v[j]->x[p], 0));
                                L->v[j]->i[p] = i;
                                p++;
                            }*/
                        }
                    }
                    printf("%s;%%   --> %ld(%ld) \n",j==An-1?"]'":"",j,p);
                    L->v[j]->nz = p;
                    Anz += p;
                }
                printf("\nu0=[\n");
                TEST_CHECK(SPEX_matrix_alloc(&U, An, An, true));
                for (j = 0; j < An; j++)
                {
                    TEST_CHECK(SPEX_vector_realloc(U->v[j], An-j));
                    TEST_CHECK(SPEX_mpz_set_si(U->v[j]->x[0], 1)); // diagnal entry
                    U->v[j]->i[0] = j;
                    p = 1;
                    for (i = 0; i < j; i++){printf("0 ");}
                            printf("1 ");
                    for (i = j+1; i < An; i++)
                    {
                        if (rand() > RAND_MAX/2)
                        {
                            printf("1 ");
                            TEST_CHECK(SPEX_mpz_set_si(U->v[j]->x[p], 1));
                            U->v[j]->i[p] = i;
                            p++;
                        }
                        else
                        {
                            printf("0 ");
                        }
                    }
                    printf("%s;%%   --> %ld(%ld) \n",j==An-1?"]":"",j,p);
                    U->v[j]->nz = p;
                    Anz += p;
                }

                k = rand()%An;
                printf("A0=l0*u0;\nk=%ld;%%1-based, for matlab code\nvk=[",k+1);

                TEST_CHECK(SPEX_vector_alloc(&vk, An, true));
                vk_nz = 0;
                for (i = 0; i < An; i++)
                {
                    if (i == k || rand() > RAND_MAX/2)
                    {
                        TEST_CHECK(SPEX_mpz_set_si(vk->x[vk_nz], 1));
                        vk->i[vk_nz] = i;
                        vk_nz++;
                        printf("1 ");
                    }
                    else
                    {
                        printf("0 ");
                    }
                }
                printf("]';\n\n\n");
                vk->nz = vk_nz;

                // printf L and U as F to file
                fprintf(mat_file, "%ld %ld %ld %ld\n", An, Anz-An, k, vk_nz);
                for (i = 0; i < An; i++)
                {
                    for (int64_t Lp = 0; Lp < L->v[i]->nz; Lp++)
                    {
                        j = L->v[i]->i[Lp];
                        SPEX_CHECK(SPEX_mpz_sgn(&sgn, L->v[i]->x[Lp]));
                        fprintf(mat_file, "%ld %ld %s\n", j,i,(sgn==0)?"0":"1");
                    }

                    for (int64_t Up = 0; Up < U->v[i]->nz; Up++)
                    {
                        j = U->v[i]->i[Up];
                        // skip the diagonal
                        if (j == i) {continue;}
                        SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[i]->x[Up]));
                        fprintf(mat_file, "%ld %ld %s\n", i,j,(sgn==0)?"0":"1");
                    }
                }

                for (j = 0; j < vk_nz; j++)
                {
                    i = vk->i[j];
                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, vk->x[j]));
                    fprintf(mat_file, "%ld %s\n", i,(sgn==0)?"0":"1");
                }
            }

            // collapse L and U
            for (i = 0; i < An; i++)
            {
                SPEX_CHECK(SPEX_vector_realloc(L->v[i], L->v[i]->nz));
                SPEX_CHECK(SPEX_vector_realloc(U->v[i], U->v[i]->nz));
            }

            // get matrix A from L*U, since D is identity matrix
            TEST_CHECK(SPEX_matrix_alloc(&A, An, An, true));
            for (i = 0; i < An; i++)
            {
                TEST_CHECK(SPEX_vector_realloc(A->v[i], An));
                for (p = 0; p < An; p++)
                {
                    TEST_CHECK(SPEX_mpz_set_si(A->v[i]->x[p], 0));
                    A->v[i]->i[p] = p;
                }
                A->v[i]->nz = An;
            }
            for (int64_t iter = 0; iter < An; iter++)
            {
                for (int64_t Up = 0; Up < U->v[iter]->nz; Up++)
                {
                    SPEX_CHECK(SPEX_mpz_sgn(&sgn, U->v[iter]->x[Up]));
                    if (sgn == 0) {continue;}
                    j = U->v[iter]->i[Up];
                    for (int64_t Lp = 0; Lp < L->v[iter]->nz; Lp++)
                    {
                        i = L->v[iter]->i[Lp];
                        // A(i,j) += L(i,iter)*U(iter,j)
                        TEST_CHECK(SPEX_mpz_addmul(A->v[j]->x[i],
                                        U->v[iter]->x[Up], L->v[iter]->x[Lp]));
                    }
                }
            }
            fclose(mat_file); mat_file = NULL;

            // create d, sd as ones(An,1), S as ones(An,3),
            // P, Q, P_inv, Q_inv as 0:An-1 
            d  = spex_create_mpz_array(An);
            sd = spex_create_mpz_array(An);
            S  = spex_create_mpq_array(3*An);
            P     = (int64_t*) SPEX_malloc(An*sizeof(int64_t));
            Q     = (int64_t*) SPEX_malloc(An*sizeof(int64_t));
            P_inv = (int64_t*) SPEX_malloc(An*sizeof(int64_t));
            Q_inv = (int64_t*) SPEX_malloc(An*sizeof(int64_t));
            if (!d || !sd || !S || !P || !Q || !P_inv || !Q_inv)
            {
                SPEX_FREE_ALL;
                return 0;
            }
            for (i = 0; i < An; i++)
            {
                TEST_CHECK(SPEX_mpz_set_ui(d[i],  1));
                TEST_CHECK(SPEX_mpz_set_ui(sd[i], 1));
                TEST_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 1, i), 1, 1));
                TEST_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 2, i), 1, 1));
                TEST_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 3, i), 1, 1));
                P[i] = i;
                Q[i] = i;
                P_inv[i] = i;
                Q_inv[i] = i;
            }
            info=SPEX_LUU(A, L, U, d, sd, S, P, P_inv, Q, Q_inv, &vk, k, NULL);
            SPEX_PRINT_INFO(info);
            if(info!=SPEX_SINGULAR && info!=SPEX_OK){TEST_CHECK(info);}
            if (info == SPEX_SINGULAR)
            {
                printf("inserting vk would cause singularity!!!\n");
            }
            GOTCHA;

            /*TEST_CHECK(SPEX_mpz_init(tmpz));
            TEST_CHECK(SPEX_mpq_init(tmpq));
            for (i = 0; i < An; i++)
            {
                TEST_CHECK(SPEX_mpq_mul(tmpq, SPEX_2D(S, 1, i), SPEX_2D(S, 3, i)));

                for (p = 0; p < L->v[i]->nz; p++)
                {
                    j = L->v[i]->i[p];
                    TEST_CHECK(SPEX_mpz_divexact(tmpz, L->v[i]->x[p],
                        SPEX_MPQ_DEN(tmpq)));
                    TEST_CHECK(SPEX_mpz_mul(tmpz, tmpz, SPEX_MPQ_NUM(tmpq)));
                    TEST_CHECK(SPEX_gmp_printf("(%ld)%Zd ",j,tmpz));
                }
                TEST_CHECK(SPEX_gmp_printf("......col %ld(S=%Qd)\n",i,tmpq));
            }
            for (i = 0; i < An; i++)
            {
                TEST_CHECK(SPEX_mpq_mul(tmpq, SPEX_2D(S, 2, i), SPEX_2D(S, 3, i)));

                for (p = 0; p < U->v[i]->nz; p++)
                {
                    j = U->v[i]->i[p];
                    TEST_CHECK(SPEX_mpz_divexact(tmpz, U->v[i]->x[p],
                        SPEX_MPQ_DEN(tmpq)));
                    TEST_CHECK(SPEX_mpz_mul(tmpz, tmpz, SPEX_MPQ_NUM(tmpq)));
                    TEST_CHECK(SPEX_gmp_printf("(%ld)%Zd ",j,tmpz));
                }
                TEST_CHECK(SPEX_gmp_printf("......col %ld(S=%Qd)\n",i,tmpq));
            }
            printf("P=[%ld %ld %ld %ld]\n",P[0], P[1], P[2], P[3]);
            printf("Q=[%ld %ld %ld %ld]\n",Q[0], Q[1], Q[2], Q[3]);*/
            SPEX_FREE_ALL;
            GOTCHA;
        }
    }

    return 0;
}

