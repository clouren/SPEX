//------------------------------------------------------------------------------
// Tcov/tcov_for_update.c: test coverage for SPEX_Update
//------------------------------------------------------------------------------

// SPEX: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX/License for the license.

//------------------------------------------------------------------------------

/*
 * test coverage: For each iteration in this test, a n-by-n sparse matrices L
 *                and U with nonzero entries be all 1 are generated, and
 *                A=LD^(-1)U is computed correspondingly. Then a n-by-1 sparse
 *                vector to be swapped with the k-th(randomly selected) column
 *                of A is randomly generated. Finally, all these components are
 *                used as input for the SPEX_Update_LU_ColRep for test. For
 *                each given n, the above process will run for 100 times
 *
 * usage:
 * *****************test for LU update for column replacement*******************
 * 1) ./tcov_for_update 0
 *                This will have n default as 4. Additional inner loop will
 *                iterate for malloc_count initialized from 0 to
 *                MAX_MALLOC_COUNT, break when malloc_count>0 at the end of
 *                inner loop.
 * 2) ./tcov_for_update 0 n
 *                This will have n set as input. n must be greater than 0.
 *                Otherwise, it will be the usages 4).  Additional inner loop
 *                will iterate for malloc_count initialized from 0 to
 *                MAX_MALLOC_COUNT, break when malloc_count>0 at the end of
 *                inner loop.
 * 3) ./tcov_for_update 0 n_min n_max
 *                This will have n iterate from n_min to n_max-1. n must be
 *                greater than 0. Otherwise, it will be the usage 5).
 *                Additional inner loop will iterate for malloc_count
 *                initialized from 0 to MAX_MALLOC_COUNT, break when
 *                malloc_count>0 at the end of inner loop.
 * 4) ./tcov_for_update 0 0
 *                This will use Mats4Tcov/SPEX_Update/mat1.txt as input to
 *                obtain L, U, vk and k to perform the test for 1 iteration.
 * 5) ./tcov_for_update 0 0 max_file_index
 *                This will use mat1.txt .... mat[max_file_index].txt in
 *                Mats4Tcov/SPEX_Update as input to obtain L, U, vk, and k.
 *                e.g., ./tcov_for_update 0 3 will use mat1.txt mat2.txt and
 *                mat3.txt as input in each iteration.
 * 6) ./tcov_for_update 0 0 max_file_num file1 file2 ... file[max_file_num]
 *                This will use the given files (total number is
 *                max_file_number and each file name is correspondingly
 *                specified as file1 etc) as input to obtain L, U, vk and k to
 *                perform the test.
 *
 * ****************test for Cholesky Rank-1 Update/Downdate*********************
 * 7) ./tcov_for_update 1
 *
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
    if (!pretend_to_fail && mat_file != NULL)    \
    {fclose(mat_file); mat_file = NULL;}         \
    TEST_OK(SPEX_matrix_free(&L, option));       \
    TEST_OK(SPEX_matrix_free(&U, option));       \
    TEST_OK(SPEX_matrix_free(&A, option));       \
    TEST_OK(SPEX_matrix_free(&b, option));       \
    TEST_OK(SPEX_matrix_free(&tmpA, option));    \
    TEST_OK(SPEX_matrix_free(&vk, option));      \
    TEST_OK(SPEX_matrix_free(&rhos, option));    \
    SPEX_FREE(h);                                \
    SPEX_FREE(P);                                \
    SPEX_FREE(P_inv);                            \
    SPEX_FREE(Q);                                \
    SPEX_FREE(Q_inv);                            \
    if (F) SPEX_MPQ_CLEAR(F->scale_for_A);       \
    SPEX_FREE (F) ;                              \
    SPEX_MPZ_CLEAR(tmpz);                        \
    SPEX_MPQ_CLEAR(tmpq);                        \
    SPEX_FREE(option);                           \
    TEST_OK(SPEX_finalize()) ;                   \
}

#define MY_PR(...) { if (pr >= 1) printf(__VA_ARGS__); }

#include "tcov_malloc_test.h"
#include "simple_rand.h"


#define MAX_MALLOC_COUNT 1000000

int64_t Aold[16] = {3, 5, 6, 7, 11, 2, -7, 10, 8, 3, -2, -2, 7, 5, 1, -6};//CSC
int64_t Lold[10] = {3, 5, 6, 7, -49, -87, -47, -17, 527, 884};// CSC
int64_t Uold[10] = {3, 11, 8, 7, -49, -31, -20, -17, 57,884};//CSR
int64_t Ak_new[4]= {1,4,7,11}; // new column
#include <assert.h>

int main( int argc, char* argv[])
{
    SPEX_info info;
    int sgn;
    int test_type = 0; //0: Update LU column replacement; 1: rank-1 Chol
    int pr = 1; // 1: print; 0: don't print
    bool pretend_to_fail = false;
    int64_t n, n_min = 4, n_max, An=0, Anz, vk_nz, i, j, k = 0, p;
    int64_t max_outer_loop = 1, max_inner_loop = 1;
    uint64_t seed = 1; // seed for random number generator
    bool read_matrix = true;

    if (argc >= 2)
    {
        test_type = atoi(argv[1]);
        if (test_type < 0 || test_type > 1)
        {
            printf("test type must be in [0,1]\n");
            return 0;
        }
    }

    if (argc >= 3)
    {
        n_min = atoi(argv[2]);
        if (n_min < 0) {return 0;}
    }
    n_max = n_min+1;
    if (argc >= 4)
    {
        n_max = atoi(argv[3]);
    }
    if (n_min == 0)
    {
        read_matrix = true;
        max_outer_loop =1;
        max_inner_loop = 1;
        if (argc > 4 && argc < n_max+4)
        {
            printf("not enough file provided! require number:%ld\n", n_max);
            return 0;
        }
    }
    else
    {
        read_matrix = false;
        max_outer_loop = 10;
        max_inner_loop = 10;
    }

    for (n = n_min; n < n_max; n++)
    {
        for (int64_t loop = 0; loop < max_outer_loop ; loop++)
        {
            for (int64_t kk = 0; kk < MAX_MALLOC_COUNT; kk++)
            {
                if (read_matrix)
                {
                    kk = MAX_MALLOC_COUNT;
                }
                malloc_count = kk;
                pretend_to_fail = false;
                //--------------------------------------------------------------
                // Initialize SPEX Update process
                //--------------------------------------------------------------

                TEST_OK(SPEX_initialize_expert (tcov_malloc, tcov_calloc,
                    tcov_realloc, tcov_free)) ;
                if (pretend_to_fail) {continue;}

                //--------------------------------------------------------------
                // Allocate memory
                //--------------------------------------------------------------

                SPEX_options* option = NULL;
                SPEX_factorization *F = NULL;
                SPEX_matrix *A = NULL, *tmpA = NULL, *vk = NULL, *b = NULL;
                SPEX_matrix *L = NULL, *U = NULL, *rhos = NULL;
                int64_t *P = NULL, *P_inv = NULL, *Q = NULL, *Q_inv = NULL;
                int64_t *h = NULL;
                mpz_t tmpz; SPEX_MPZ_SET_NULL(tmpz);
                mpq_t tmpq; SPEX_MPQ_SET_NULL(tmpq);
                FILE *mat_file = NULL;
                TEST_CHECK(SPEX_create_default_options(&option));
                if (pretend_to_fail) {continue;}
                TEST_CHECK(SPEX_mpz_init(tmpz));
                if (pretend_to_fail) {continue;}
                TEST_CHECK(SPEX_mpq_init(tmpq));
                if (pretend_to_fail) {continue;}
                simple_rand_seed(seed);
                option->print_level = 3;

                // allocate memory space for the factorization
                F = (SPEX_factorization*) SPEX_calloc(1,
                    sizeof(SPEX_factorization));
                if (F == NULL) {TEST_CHECK(SPEX_OUT_OF_MEMORY); continue;}
                // set factorization kind
                F->kind = (test_type == 0) ?
                    SPEX_LU_FACTORIZATION : SPEX_CHOLESKY_FACTORIZATION;
                // set F as updatable
                F->updatable = true;
                // Allocate and set scale_for_A
                TEST_CHECK(SPEX_mpq_init(F->scale_for_A));
                if (pretend_to_fail) {continue;}
                TEST_CHECK(SPEX_mpq_set_ui (F->scale_for_A, 1, 1));
                if (pretend_to_fail) {continue;}

                if (read_matrix)
                {
                    // read in L and U
                    if (argc > 4)
                    {
                        mat_file = fopen(argv[4+n], "r");
                    }
                    else
                    {
                        char mat_name[100];
                        sprintf(mat_name,"Mats4Tcov/SPEX_Update/mat%ld.txt",
                            n+1);
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

                    // allocate space for L and U
                    TEST_CHECK(SPEX_matrix_allocate(&(L), SPEX_DYNAMIC_CSC,
                        SPEX_MPZ, An, An, 0, false, true, option));
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK(SPEX_matrix_allocate(&(U), SPEX_DYNAMIC_CSC,
                        SPEX_MPZ, An, An, 0, false, true, option));
                    if (pretend_to_fail) {continue;}
                    for (j = 0; j < An && !pretend_to_fail; j++)
                    {
                        TEST_CHECK(SPEX_vector_realloc(L->v[j], An-j, option));
                        if (pretend_to_fail) {break;}
                        TEST_CHECK(SPEX_vector_realloc(U->v[j], An-j, option));
                        if (pretend_to_fail) {break;}
                    }
                    if (pretend_to_fail) {continue;}

                    // allocate space for vk
                    TEST_CHECK(SPEX_matrix_allocate(&vk, SPEX_DYNAMIC_CSC,
                        SPEX_MPZ, An, 1, 0, false, true, option));
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK(SPEX_vector_realloc(vk->v[0], An, option));
                    if (pretend_to_fail) {continue;}

                    // Read in the values from file
                    int64_t tmp;
                    for (int64_t ai = 0; ai < Anz && !pretend_to_fail; ai++)
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
                            if (pretend_to_fail) {break;}
                            L->v[j]->i[L->v[j]->nz] = i;
                            L->v[j]->nz++;
                        }
                        if (j >= i)
                        {
                            TEST_CHECK(SPEX_mpz_set_si(
                                U->v[i]->x[U->v[i]->nz], tmp));
                            if (pretend_to_fail) {break;}
                            U->v[i]->i[U->v[i]->nz] = j;
                            U->v[i]->nz++;
                        }
                    }
                    if (pretend_to_fail) {continue;}

                    for (j = 0; j < vk_nz && !pretend_to_fail; j++)
                    {
                        s = fscanf(mat_file, "%"PRId64" %"PRId64"\n", &i, &tmp);
                        if ((feof(mat_file) && j != vk_nz-1) || s < 2)
                        {
                            printf ("premature end-of-file\n") ;
                            SPEX_FREE_ALL;
                            return 0;
                        }
                        TEST_CHECK(SPEX_mpz_set_si(vk->v[0]->x[j], tmp));
                        if (pretend_to_fail) {break;}
                        vk->v[0]->i[j] = i;
                    }
                    if (pretend_to_fail) {continue;}
                    vk->v[0]->nz = vk_nz;
                }
                else // randomly generate L, U, vk and k
                {
                    MY_PR("\n\n\n+++++++++++++++++++++++++++++++++++++++++\n");
                    MY_PR("++++++++++++++case %ld++++++++++++++++++++\n",loop);
                    MY_PR("+++++++++++++++++++++++++++++++++++++++++++\n\n\n");
                    An = n;
                    Anz = 0;
                    char *mat_name = "Mats4Tcov/SPEX_Update/mat.txt";
                    mat_file = fopen(mat_name, "w");
                    if (mat_file == NULL)
                    {
                        perror("Error while opening the file");
                        SPEX_FREE_ALL;
                        return 0;
                    }

                    // randomly generate L, U and vk
                    TEST_CHECK(SPEX_matrix_allocate(&L, SPEX_DYNAMIC_CSC,
                        SPEX_MPZ, An, An, 0, false, true, option));
                    if (pretend_to_fail) {continue;}

                    MY_PR("l0=[\n");
                    for (j = 0; j < An && !pretend_to_fail; j++)
                    {
                        TEST_CHECK(SPEX_vector_realloc(L->v[j], An-j, option));
                        if (pretend_to_fail) {break;}
                        // diagnal entry
                        TEST_CHECK(SPEX_mpz_set_si(L->v[j]->x[0], 1));
                        if (pretend_to_fail) {break;}
                        L->v[j]->i[0] = j;
                        p = 1;
                        for (i = 0; i < j; i++){MY_PR("0 ");}
                                MY_PR("1 ");
                        for (i = j+1; i < An && !pretend_to_fail; i++)
                        {
                            if (simple_rand() > SIMPLE_RAND_MAX/2)
                            {
                                MY_PR("1 ");
                                TEST_CHECK(SPEX_mpz_set_si(L->v[j]->x[p], 1));
                                if (pretend_to_fail) {break;}
                                L->v[j]->i[p] = i;
                                p++;
                            }
                            else
                            {
                                MY_PR("0 ");
                                /*if (simple_rand() > SIMPLE_RAND_MAX/10)
                                {
                                    TEST_CHECK(SPEX_mpz_set_si(
                                        L->v[j]->x[p], 0));
                                    if (pretend_to_fail) {break;}
                                    L->v[j]->i[p] = i;
                                    p++;
                                }*/
                            }
                        }
                        if (pretend_to_fail) {break;}
                        MY_PR("%s;%%   --> %ld(%ld) \n",j==An-1?"]'":"",j,p);
                        L->v[j]->nz = p;
                        Anz += p;
                    }
                    if (pretend_to_fail) {continue;}
                    MY_PR("\nu0=[\n");
                    TEST_CHECK(SPEX_matrix_allocate(&U, SPEX_DYNAMIC_CSC,
                        SPEX_MPZ, An, An, 0, false, true, option));
                    if (pretend_to_fail) {continue;}
                    for (j = 0; j < An && !pretend_to_fail; j++)
                    {
                        TEST_CHECK(SPEX_vector_realloc(U->v[j], An-j, option));
                        if (pretend_to_fail) {break;}
                        // diagnal entry
                        TEST_CHECK(SPEX_mpz_set_si(U->v[j]->x[0], 1));
                        if (pretend_to_fail) {break;}
                        U->v[j]->i[0] = j;
                        p = 1;
                        for (i = 0; i < j; i++){MY_PR("0 ");}
                                MY_PR("1 ");
                        for (i = j+1; i < An && !pretend_to_fail; i++)
                        {
                            if (simple_rand() > SIMPLE_RAND_MAX/2)
                            {
                                MY_PR("1 ");
                                TEST_CHECK(SPEX_mpz_set_si(U->v[j]->x[p], 1));
                                if (pretend_to_fail) {break;}
                                U->v[j]->i[p] = i;
                                p++;
                            }
                            else
                            {
                                MY_PR("0 ");
                            }
                        }
                        if (pretend_to_fail) {break;}
                        MY_PR("%s;%%   --> %ld(%ld) \n",j==An-1?"]":"",j,p);
                        U->v[j]->nz = p;
                        Anz += p;
                    }
                    if (pretend_to_fail) {continue;}

                    k = simple_rand()%An;
                    MY_PR("A0=l0*u0;\nk=%ld;%%1-based, for matlab code\nvk=[",
                        k+1);

                    // allocate space for vk
                    TEST_CHECK(SPEX_matrix_allocate(&vk, SPEX_DYNAMIC_CSC,
                        SPEX_MPZ, An, 1, 0, false, true, option));
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK(SPEX_vector_realloc(vk->v[0], An, option));
                    if (pretend_to_fail) {continue;}
                    vk_nz = 0;
                    for (i = 0; i < An && !pretend_to_fail; i++)
                    {
                        if (i == k || simple_rand() > SIMPLE_RAND_MAX/2)
                        {
                            TEST_CHECK(SPEX_mpz_set_si(vk->v[0]->x[vk_nz], 1));
                            if (pretend_to_fail) {break;}
                            vk->v[0]->i[vk_nz] = i;
                            vk_nz++;
                            MY_PR("1 ");
                        }
                        else
                        {
                            MY_PR("0 ");
                        }
                    }
                    if (pretend_to_fail) {continue;}
                    MY_PR("]';\n\n\n");
                    vk->v[0]->nz = vk_nz;

                    // print L and U as F to file
                    fprintf(mat_file, "%ld %ld %ld %ld\n", An, Anz-An, k,vk_nz);
                    for (i = 0; i < An && !pretend_to_fail; i++)
                    {
                        for (int64_t Lp = 0; Lp < L->v[i]->nz &&
                            !pretend_to_fail; Lp++)
                        {
                            j = L->v[i]->i[Lp];
                            TEST_CHECK(SPEX_mpz_sgn(&sgn, L->v[i]->x[Lp]));
                            if (pretend_to_fail) {break;}
                            fprintf(mat_file, "%ld %ld %s\n", j, i,
                                (sgn==0)?"0":"1");
                        }
                        if (pretend_to_fail) {break;}

                        for (int64_t Up = 0; Up < U->v[i]->nz &&
                            !pretend_to_fail; Up++)
                        {
                            j = U->v[i]->i[Up];
                            // skip the diagonal
                            if (j == i) {continue;}
                            TEST_CHECK(SPEX_mpz_sgn(&sgn, U->v[i]->x[Up]));
                            if (pretend_to_fail) {break;}
                            fprintf(mat_file, "%ld %ld %s\n", i, j,
                                (sgn==0)?"0":"1");
                        }
                        if (pretend_to_fail) {break;}
                    }
                    if (pretend_to_fail) {continue;}

                    for (j = 0; j < vk_nz && !pretend_to_fail; j++)
                    {
                        i = vk->v[0]->i[j];
                        TEST_CHECK(SPEX_mpz_sgn(&sgn, vk->v[0]->x[j]));
                        if (pretend_to_fail) {break;}
                        fprintf(mat_file, "%ld %s\n", i,(sgn==0)?"0":"1");
                    }
                    if (pretend_to_fail) {continue;}
                }

                // collapse L and U
                for (i = 0; i < An && !pretend_to_fail; i++)
                {
                    TEST_CHECK(SPEX_vector_realloc(L->v[i], L->v[i]->nz,
                        option));
                    if (pretend_to_fail) {break;}
                    TEST_CHECK(SPEX_vector_realloc(U->v[i], U->v[i]->nz,
                        option));
                    if (pretend_to_fail) {break;}
                }
                if (pretend_to_fail) {continue;}

                // get matrix A from L*U, since D is identity matrix
                TEST_CHECK(SPEX_matrix_allocate(&A, SPEX_DYNAMIC_CSC, SPEX_MPZ,
                    An, An, 0, false, true, option));
                if (pretend_to_fail) {continue;}
                for (i = 0; i < An && !pretend_to_fail; i++)
                {
                    TEST_CHECK(SPEX_vector_realloc(A->v[i], An, option));
                    if (pretend_to_fail) {break;}
                    for (p = 0; p < An && !pretend_to_fail; p++)
                    {
                        TEST_CHECK(SPEX_mpz_set_si(A->v[i]->x[p], 0));
                        if (pretend_to_fail) {break;}
                        A->v[i]->i[p] = p;
                    }
                    if (pretend_to_fail) {break;}
                    A->v[i]->nz = An;
                }
                if (pretend_to_fail) {continue;}

                for (int64_t iter = 0; iter < An && !pretend_to_fail; iter++)
                {
                    for (int64_t Up = 0; Up < U->v[iter]->nz &&
                        !pretend_to_fail; Up++)
                    {
                        TEST_CHECK(SPEX_mpz_sgn(&sgn, U->v[iter]->x[Up]));
                        if (pretend_to_fail) {break;}
                        if (sgn == 0) {continue;}
                        j = U->v[iter]->i[Up];
                        for (int64_t Lp = 0; Lp < L->v[iter]->nz &&
                            !pretend_to_fail; Lp++)
                        {
                            i = L->v[iter]->i[Lp];
                            // A(i,j) += L(i,iter)*U(iter,j)
                            TEST_CHECK(SPEX_mpz_addmul(A->v[j]->x[i],
                                       U->v[iter]->x[Up], L->v[iter]->x[Lp]));
                            if (pretend_to_fail) {break;}
                        }
                        if (pretend_to_fail) {break;}
                    }
                    if (pretend_to_fail) {break;}
                }
                if (pretend_to_fail) {continue;}
                fclose(mat_file); mat_file = NULL;

                // create rhos as ones(An,1), and P, Q, P_inv, Q_inv as 0:An-1 
                P     = (int64_t*) SPEX_malloc(An*sizeof(int64_t));
                if (P     == NULL) {TEST_CHECK(SPEX_OUT_OF_MEMORY); continue;}
                Q     = (int64_t*) SPEX_malloc(An*sizeof(int64_t));
                if (Q     == NULL) {TEST_CHECK(SPEX_OUT_OF_MEMORY); continue;}
                P_inv = (int64_t*) SPEX_malloc(An*sizeof(int64_t));
                if (P_inv == NULL) {TEST_CHECK(SPEX_OUT_OF_MEMORY); continue;}
                Q_inv = (int64_t*) SPEX_malloc(An*sizeof(int64_t));
                if (Q_inv == NULL) {TEST_CHECK(SPEX_OUT_OF_MEMORY); continue;}

                TEST_CHECK(SPEX_matrix_allocate(&rhos, SPEX_DENSE, SPEX_MPZ,
                    An, 1, An, false, true, option));
                if (pretend_to_fail) {continue;}
                for (i = 0; i < An && !pretend_to_fail; i++)
                {
                    TEST_CHECK(SPEX_mpz_set_ui(rhos->x.mpz[i], 1));
                    if (pretend_to_fail) {break;}
                    P[i] = i;
                    Q[i] = i;
                    P_inv[i] = i;
                    Q_inv[i] = i;
                }
                if (pretend_to_fail) {continue;}

                // perform LU update for column replacement
                bool got_singular_matrix = false;
                for (int64_t inner = 0; inner < max_inner_loop &&
                    !pretend_to_fail; inner++)
                {
                    if (test_type == 0)
                    {
                        F->L = L;
                        F->U = U;
                        F->rhos = rhos;
                        F->P_perm = P;
                        F->Pinv_perm = P_inv;
                        F->Q_perm = Q;
                        F->Qinv_perm = Q_inv;
                        info = SPEX_Update_LU_ColRep(F, vk, k, option);
                        //info = SPEX_Update_LU_ColRep(A, L, U, rhos, P, P_inv, Q,
                          //  Q_inv, &vk, k, option);
                    }
                    else if (test_type == 1)
                    {
                        F->L = L;
                        F->rhos = rhos;
                        F->P_perm = P;
                        F->Pinv_perm = P_inv;
                        // always perform update in the first loop
                        int64_t sigma = 1;
                        if (inner > 0)
                        {
                            sigma = (simple_rand() > SIMPLE_RAND_MAX/2) ? -1: 1;
                        }
                        printf("performing Chol rank-1 %s\n", sigma < 0 ?
                            "downdate" : "update");
                        // always shrink vk
                        // TODO delete?
                        /*TEST_CHECK(SPEX_vector_realloc(vk, vk->nz, option));
                        if (pretend_to_fail) {break;}*/
                        info = SPEX_Update_Chol_Rank1(F, vk, sigma, option);
                        //(L, rhos, P, P_inv, vk, sigma, option);
                    }
                        
                    if (info == SPEX_SINGULAR)
                    {
                        printf("This update would cause singularity!!!\n");
                        got_singular_matrix = true;
                        break;
                    }
                    else if (info != SPEX_SINGULAR && info != SPEX_OK)
                    {
                        TEST_CHECK(info);
                        if (pretend_to_fail) {break;}
                    }

                    // randomly generate a new vk and perform another update
                    if (read_matrix)       { break; }

                    k = simple_rand()%An;
                    MY_PR("k=%ld;%%1-based, for matlab code\nvk=[",
                        k+1);

                    // TODO delete?
                    /*if (test_type == 1)
                    {
                        // reallocate An space for vk
                        TEST_CHECK(SPEX_vector_realloc(vk, An, option));
                        if (pretend_to_fail) {break;}
                    }
                    vk_nz = 0;
                    for (i = 0; i < An && !pretend_to_fail; i++)
                    {
                        if (i == k || simple_rand() > SIMPLE_RAND_MAX/2)
                        {
                            TEST_CHECK(SPEX_mpz_set_si(vk->x[vk_nz], 1));
                            if (pretend_to_fail) {break;}
                            vk->i[vk_nz] = i;
                            vk_nz++;
                            MY_PR("1 ");
                        }
                        else
                        {
                            MY_PR("0 ");
                        }
                    }
                    if (pretend_to_fail) {break;}
                    MY_PR("]';\n\n\n");
                    vk->nz = vk_nz;*/
                }
                if (pretend_to_fail) {continue;}
                if (got_singular_matrix)
                {
                    SPEX_FREE_ALL;
                    if (malloc_count >= 0)
                    {
                        printf("require malloc count: %ld\n\n\n\n",
                            kk-malloc_count);
                        seed = simple_rand_getseed();
                        break;
                    }
                }

                // FIXME
                //print U which will apply scale in each vector
                /*TEST_CHECK(SPEX_matrix_check(F->U, option));
                if (pretend_to_fail) {continue;}
                // get the static F
                TEST_CHECK(SPEX_factorization_convert(F, option));
                if (pretend_to_fail) {continue;}*/

                if (read_matrix)
                {
                    //----------------------------------------------------------
                    // test some other functions
                    //----------------------------------------------------------
                    // test SPEX_matrix_copy and SPEX_transpose
                    for (int kind = 0; kind < 4 && !pretend_to_fail; kind++)
                    {
                        int max_type = 5;
                        if (kind == 3)
                        {
                            max_type = 1; // SPEX_MPZ
                        }
                        for (int type = 0; type < max_type && !pretend_to_fail;
                             type++)
                        {
                            // convert A to kind x type matrix tmpA
                            TEST_OK(SPEX_matrix_free(&tmpA, option));
                            if (pretend_to_fail) {break;}
                            TEST_CHECK(SPEX_matrix_copy(&tmpA, (SPEX_kind) kind,
                                (SPEX_type) type, A, option));
                            if (pretend_to_fail) {break;}

                            if (kind == 0)
                            {
                                // if kind == SPEX_CSC, transpose the resulted
                                // matrix
                                TEST_OK(SPEX_matrix_free(&A, option));
                                if (pretend_to_fail) {break;}
                                TEST_CHECK(SPEX_transpose(&A, tmpA, option));
                                if (pretend_to_fail) {break;}
                            }

                            // convert tmpA back to original A
                            TEST_OK(SPEX_matrix_free(&A, option));
                            if (pretend_to_fail) {break;}
                            TEST_CHECK(SPEX_matrix_copy(&A, SPEX_DYNAMIC_CSC,
                                SPEX_MPZ, tmpA, option));
                            if (pretend_to_fail) {break;}
                        }
                        if (pretend_to_fail) {break;}

                    }
                    if (pretend_to_fail) {continue;}

                    //----------------------------------------------------------
                    // failure cases
                    //----------------------------------------------------------
                    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                    // fail SPEX_Update_Solve
                    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                    TEST_CHECK_FAILURE(SPEX_Update_solve(NULL, NULL, NULL,
                        option), SPEX_INCORRECT_INPUT);
                    if (pretend_to_fail) {continue;}

#if 0
                    // intentionally incorrect use of matrix_canonicalize
                    // to produce a bad L so as to fail spex_update_ipge, which
                    // is called by SPEX_Update_Solve
                    TEST_CHECK(SPEX_Update_matrix_canonicalize(L, Q_inv,
                        option));
                    if (pretend_to_fail) {continue;}
                    // make the permutation broken
                    int64_t P0 = P[0]; P[0] = P[1]; P[1] = P0;

                    // use SPEX_Update_Solve
                    h = (int64_t*) SPEX_malloc(An*sizeof(int64_t));
                    if (h == NULL)
                    {
                        TEST_CHECK(SPEX_OUT_OF_MEMORY);
                        continue;
                    }
                    TEST_CHECK(SPEX_matrix_allocate(&b, SPEX_DENSE, SPEX_MPZ,
                        An, 1, An, false, true, option));
                    if (pretend_to_fail) {continue;}
                    for (i = 0; i < An && !pretend_to_fail; i++)
                    {
                        TEST_CHECK(SPEX_mpz_set_si(b->x.mpz[i], i+1));
                        if (pretend_to_fail) {break;}
                    }
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK_FAILURE(SPEX_Update_solve(&b_sol, F, b, option),
                        SPEX_INCORRECT_INPUT);
                    if (pretend_to_fail) {continue;}
#endif

                    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                    // test and fail SPEX_matrix_check
                    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                    TEST_CHECK(SPEX_mpq_set_ui(L->v[An-1]->scale, 1, 100000));
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK_FAILURE(SPEX_matrix_check(L, option),
                        SPEX_INCORRECT_INPUT);
                    if (pretend_to_fail) {continue;}

                    // duplicate index
                    if (L->v[0]->nz > 1)
                    {
                        L->v[0]->i[1] = L->v[0]->i[0];
                        TEST_CHECK_FAILURE(SPEX_matrix_check(L, option),
                            SPEX_INCORRECT_INPUT);
                        if (pretend_to_fail) {continue;}
                    }

                    // index out of range
                    L->v[0]->i[0] = An;
                    TEST_CHECK_FAILURE(SPEX_matrix_check(L, option),
                        SPEX_INCORRECT_INPUT);
                    if (pretend_to_fail) {continue;}

                    // invalid i
                    SPEX_FREE(L->v[0]->i);
                    printf("L->v->i=%p\n",L->v[0]->i);
                    TEST_CHECK_FAILURE(SPEX_matrix_check(L, option),
                        SPEX_INCORRECT_INPUT);
                    if (pretend_to_fail) {continue;}

                    // invalid v
                    SPEX_vector **tmp_Lv = L->v;
                    L->v = NULL;
                    TEST_CHECK_FAILURE(SPEX_matrix_check(L, option),
                        SPEX_INCORRECT_INPUT);
                    if (pretend_to_fail) {continue;}
                    L->v = tmp_Lv;

                    // invalid v
                    TEST_OK(SPEX_vector_free(&(L->v[An-1]), option));
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK_FAILURE(SPEX_matrix_check(L, option),
                        SPEX_INCORRECT_INPUT);
                    if (pretend_to_fail) {continue;}

                    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                    // fail spex_update_verify
                    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                    // free L, U and A and reconstruct them as 1-by-1 matrix
                    TEST_OK(SPEX_matrix_free(&L, option));
                    if (pretend_to_fail) {continue;}
                    TEST_OK(SPEX_matrix_free(&U, option));
                    if (pretend_to_fail) {continue;}
                    TEST_OK(SPEX_matrix_free(&A, option));
                    if (pretend_to_fail) {continue;}
                    TEST_OK(SPEX_matrix_free(&rhos, option));
                    if (pretend_to_fail) {continue;}

                    // U = one(1,1)
                    TEST_CHECK(SPEX_matrix_allocate(&U, SPEX_DYNAMIC_CSC,
                        SPEX_MPZ, 1, 1, 0, false, true, option));
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK(SPEX_vector_realloc(U->v[0], 1, option));
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK(SPEX_mpz_set_ui(U->v[0]->x[0], 1));
                    if (pretend_to_fail) {continue;}
                    U->v[0]->i[0] = 0;
                    U->v[0]->nz = 1;

                    // L = one(1,1)
                    TEST_CHECK(SPEX_matrix_allocate(&L, SPEX_DYNAMIC_CSC,
                        SPEX_MPZ, 1, 1, 0, false, true, option));
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK(SPEX_vector_realloc(L->v[0], 1, option));
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK(SPEX_mpz_set_ui(L->v[0]->x[0], 1));
                    if (pretend_to_fail) {continue;}
                    L->v[0]->i[0] = 0;
                    L->v[0]->nz = 1;

                    // A = one(1,1)*2
                    TEST_CHECK(SPEX_matrix_allocate(&A, SPEX_DYNAMIC_CSC,
                        SPEX_MPZ, 1, 1, 0, false, true, option));
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK(SPEX_vector_realloc(A->v[0], 1, option));
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK(SPEX_mpz_set_ui(A->v[0]->x[0], 2));
                    if (pretend_to_fail) {continue;}
                    A->v[0]->i[0] = 0;
                    A->v[0]->nz = 1;
                    GOTCHA;

                    // rhos = ones (1,1)
                    TEST_CHECK(SPEX_matrix_allocate(&rhos, SPEX_DENSE,
                        SPEX_MPZ, 1, 1, 1, false, true, option));
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK(SPEX_mpz_set_ui(rhos->x.mpz[0], 1));
                    if (pretend_to_fail) {continue;}
                    GOTCHA;

                    // fake P and Q_inv
                    P[0] = 0;
                    Q_inv[0] = 0;
                    // FIXME move
                    F->L = L; F->U = U; F->rhos = rhos;
                    F->P_perm = P; F->Qinv_perm = Q_inv;
                    TEST_CHECK_FAILURE(spex_update_verify(F, A, option),
                        SPEX_INCORRECT);
                    if (pretend_to_fail) {continue;}
                    GOTCHA;

                    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                    // fail SPEX_vector_allocate
                    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                    SPEX_vector *tmpv = NULL;
                    TEST_CHECK_FAILURE(SPEX_vector_allocate(NULL, 1, option),
                        SPEX_INCORRECT_INPUT);
                    if (pretend_to_fail) {continue;}
                    TEST_CHECK_FAILURE(SPEX_vector_allocate(&tmpv, -1, option),
                        SPEX_INCORRECT_INPUT);
                    if (pretend_to_fail) {continue;}

                    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                    // fail SPEX_Update_LU_ColRep
                    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                    TEST_CHECK_FAILURE(SPEX_Update_LU_ColRep(NULL, NULL, k,
                        option), SPEX_INCORRECT_INPUT);
                    if (pretend_to_fail) {continue;}

                    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                    // fail SPEX_Update_Chol_Rank1
                    // . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
                    TEST_CHECK_FAILURE(SPEX_Update_Chol_Rank1(
                        NULL, NULL, 1, option), SPEX_INCORRECT_INPUT);
                    if (pretend_to_fail) {continue;}
                }

                SPEX_FREE_ALL;
                if (malloc_count >= 0)
                {
                    printf("require malloc count: %ld\n\n\n\n",kk-malloc_count);
                    seed = simple_rand_getseed();
                    break;
                }
            }
        }
    }

    printf("test finished\n");
    return 0;
}

