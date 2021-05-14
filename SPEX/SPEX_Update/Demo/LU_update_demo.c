//------------------------------------------------------------------------------
// SPEX_CHOLMOD/Demo/LU_update_demo.c: demo for SPEX_CHOLMOD library
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/*
 * test coverage
 */

#define FREE_WORKSPACE                            \
{                                                \
    if (mat_file!=NULL) fclose(mat_file);        \
    SPEX_FREE(option);                           \
    SPEX_finalize() ;                            \
}

//#include "SPEX_CHOLMOD.h"
#include "tcov_malloc_test.h"

#define OK(method)                       \
{                                                \
    info = (method) ;                            \
    if (info != SPEX_OK)                         \
    {                                            \
        FREE_WORKSPACE;                           \
        SPEX_PRINT_INFO(info);\
        return 0 ;/*continue; */                               \
    }                                            \
}


#include <assert.h>
#include <inttypes.h>
int64_t Q_for_T[177] = {4, 45, 38, 64, 161, 54, 1, 49, 78, 86, 125, 172, 57, 27, 167, 136, 6, 8, 60, 101, 148, 48, 174, 74, 41, 143, 24, 137, 113, 51, 18, 42, 97, 129, 146, 99, 106, 127, 71, 19, 144, 77, 72, 121, 46, 158, 94, 26, 10, 145, 159, 147, 156, 149, 30, 176, 111, 59, 39, 128, 157, 175, 104, 108, 169, 118, 16, 25, 56, 58, 63, 65, 70, 82, 88, 89, 92, 103, 109, 116, 120, 122, 139, 142, 162, 163, 166, 90, 5, 0, 11, 28, 29, 44, 83, 87, 124, 135, 138, 107, 165, 164, 55, 68, 75, 84, 105, 53, 12, 20, 21, 43, 47, 91, 93, 117, 119, 140, 126, 22, 23, 2, 14, 17, 36, 40, 61, 62, 79, 80, 81, 95, 141, 160, 31, 37, 52, 67, 155, 173, 35, 123, 3, 7, 9, 13, 15, 32, 33, 131, 133, 170, 171, 34, 50, 85, 96, 98, 100, 102, 110, 112, 114, 115, 130, 150, 151, 152, 168, 134, 73, 132, 153, 154, 69, 66, 76};

int64_t P_inv_for_T[177] = {78, 21, 36, 22, 65, 23, 12, 0, 53, 125, 89, 41, 45, 31, 7, 15, 25, 16, 5, 55, 66, 29, 17, 6, 34, 11, 27, 18, 168, 14, 98, 43, 26, 93, 2, 30, 122, 24, 58, 9, 94, 88, 4, 44, 39, 1, 3, 46, 10, 49, 126, 32, 13, 42, 38, 19, 33, 20, 51, 48, 54, 113, 8, 104, 171, 173, 144, 90, 57, 68, 82, 67, 52, 96, 107, 95, 71, 76, 139, 59, 62, 74, 70, 60, 80, 137, 72, 140, 35, 75, 47, 86, 92, 85, 73, 147, 103, 87, 69, 79, 91, 172, 56, 77, 40, 99, 101, 121, 108, 117, 145, 110, 118, 63, 114, 61, 162, 112, 109, 127, 115, 106, 81, 119, 157, 146, 124, 37, 123, 154, 120, 130, 163, 105, 170, 138, 161, 100, 176, 28, 141, 111, 83, 143, 132, 166, 175, 50, 135, 149, 159, 167, 142, 156, 116, 164, 158, 134, 150, 153, 160, 151, 84, 148, 155, 131, 97, 165, 152, 64, 102, 129, 128, 136, 174, 133, 169};

#define num_col_index  16
int64_t col_index[num_col_index] = {46, 12, 65, 40, 63, 9, 4, 0, 87, 16, 76, 14,54, 100, 2, 39};

int main( int argc, char* argv[])
{
    SPEX_info info;
    //------------------------------------------------------------------
    // Initialize SPEX CHOLMOD process
    //------------------------------------------------------------------

    SPEX_initialize () ;

    //------------------------------------------------------------------
    // Allocate memory
    //------------------------------------------------------------------

    int64_t m, n, nz, i, j;
    SPEX_options* option = NULL;

    SPEX_matrix *L = NULL, *U = NULL, *A = NULL, *T = NULL;
    SPEX_matrix *rhos = NULL;
    int64_t *P = NULL, *P_inv = NULL, *Q = NULL, *Q_inv = NULL, *mark = NULL;
    int64_t tmpi;
    /*mpz_t tmpz; SPEX_MPZ_SET_NULL(tmpz);
    OK(SPEX_mpz_init(tmpz));
    mpq_t tmpq; SPEX_MPQ_SET_NULL(tmpq);
    OK(SPEX_mpq_init(tmpq));*/

    char *mat_name = "../ExampleMats/10teams_mat.txt";
    FILE *mat_file = fopen(mat_name, "r");
    if (mat_file == NULL)
    {
        perror("Error while opening the file");
        FREE_WORKSPACE;
        return 0;
    }

    // Read in size of matrix & number of nonzeros
    int s = fscanf(mat_file, "%"PRId64" %"PRId64" %"PRId64"\n", &m, &n, &nz);
    if (feof(mat_file) || s < 3)
    {
        printf ("premature end-of-file\n") ;
        FREE_WORKSPACE;
        return 0;
    }

    if (m != n)
    {
        printf("input matrix must be square matrix\n");
        FREE_WORKSPACE;
        return 0;
    }

    OK(SPEX_create_default_options(&option));
    OK(SPEX_matrix_allocate(&T, SPEX_DYNAMIC_CSC, SPEX_MPZ, m, n, 0, false,
        true, option));
    for (j = 0; j < n; j++)
    {
        OK(SPEX_vector_realloc(T->v[j], m, option));
    }
    // Read in the values from file
    int64_t Tij;
    for (int64_t k = 0; k < nz; k++)
    {
        s = fscanf(mat_file, "%"PRId64" %"PRId64" %"PRId64"\n", &i, &j, &Tij);
        if ((feof(mat_file) && k != nz-1) || s < 3)
        {
            printf ("premature end-of-file\n") ;
            FREE_WORKSPACE;
            return 0;
        }

        // Conversion from 1 based to 0 based
        T->v[j-1]->i[T->v[j-1]->nz] = P_inv_for_T[i-1];
        //Tij = rand()%5+1;
        OK(SPEX_mpz_set_si(T->v[j-1]->x[T->v[j-1]->nz], Tij));
        T->v[j-1]->nz++;
    }

    OK(SPEX_matrix_allocate(&L, SPEX_DYNAMIC_CSC, SPEX_MPZ, m, n, 0, false,
        true, option));
    OK(SPEX_matrix_allocate(&U, SPEX_DYNAMIC_CSC, SPEX_MPZ, m, n, 0, false,
        true, option));
    OK(SPEX_matrix_allocate(&A, SPEX_DYNAMIC_CSC, SPEX_MPZ, m, n, 0, false,
        true, option));
    OK(SPEX_matrix_allocate(&rhos, SPEX_DENSE, SPEX_MPZ, n, 1,
	n, false, true, option));
    P     = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    Q     = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    P_inv = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    Q_inv = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    if (!P || !Q || !P_inv || !Q_inv)
    {
        FREE_WORKSPACE;
        return 0;
    }
    // set A=L=U=P=Q=I
    for (i = 0; i < n; i++)
    {
        OK(SPEX_vector_realloc(A->v[i], 1, option));
        OK(SPEX_vector_realloc(L->v[i], 1, option));
        OK(SPEX_vector_realloc(U->v[i], 1, option));

        OK(SPEX_mpz_set_ui(A->v[i]->x[0], 1));
        A->v[i]->i[0] = i;
        A->v[i]->nz = 1;
        OK(SPEX_mpz_set_ui(L->v[i]->x[0], 1));
        L->v[i]->i[0] = i;
        L->v[i]->nz = 1;
        OK(SPEX_mpz_set_ui(U->v[i]->x[0], 1));
        U->v[i]->i[0] = i;
        U->v[i]->nz = 1;
        OK(SPEX_mpz_set(rhos->x.mpz[i], L->v[i]->x[0]));
        P[i] = i;
        Q[i] = i;
        P_inv[i] = i;
        Q_inv[i] = i;
    }

    mark = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    for (i = 0; i < n; i++)
    {
        mark[i] = -1;
    }

    for (int64_t k = 0; k < num_col_index; k++)
    {
        j = col_index[k];
        printf("-----------------------------------------------\n");
        printf("j = %ld\n", j);
        info = SPEX_Update_LU_ColRep(A, L, U, rhos, P, P_inv, Q, Q_inv,
	    &(T->v[Q_for_T[j]]), j, option);
        if (info != SPEX_OK)
        {
            printf("k=%ld\n",k);
            return 0;
        }

        mark[j] = k;
    }

    int64_t found_index = num_col_index-1;

    // a random process to find a column from T and use that to replace
    // the corresponding column in A, then update the LU factorization of A
    while (found_index < n)
    {
        tmpi = 0;
        while (mark[tmpi] > -1 || mark[tmpi] == -found_index)
        {
            tmpi = rand()%n;
        }
        j = tmpi;
        printf("-----------------------------------------------\n");
        printf("j = %ld mark=%ld found_index = %ld\n", j, mark[j],found_index);


        info = SPEX_Update_LU_ColRep(A, L, U, rhos, P, P_inv, Q, Q_inv,
	    &(T->v[Q_for_T[j]]), j, option);
        if (info != SPEX_SINGULAR && info != SPEX_OK)
        {
            FREE_WORKSPACE;
            return 0;
        }
        else if (info == SPEX_OK)
        {
            found_index ++;
            mark[j] = found_index;
        }
        else
        {
            printf("resulting singular matrix with column %ld\n",j);
            mark[j] = -found_index;
        }

    /*    for (i = 0; i < n; i++)
        {
            OK(SPEX_mpq_mul(tmpq, SPEX_2D(S, 1, i), SPEX_2D(S, 3, i)));

            for (p = 0; p < L->v[i]->nz; p++)
            {
                j = L->v[i]->i[p];
                OK(SPEX_mpz_divexact(tmpz, L->v[i]->x[p],
                    SPEX_MPQ_DEN(tmpq)));
                OK(SPEX_mpz_mul(tmpz, tmpz, SPEX_MPQ_NUM(tmpq)));
                OK(SPEX_gmp_printf("(%ld)%Zd ",j,tmpz));
            }
            OK(SPEX_gmp_printf("......col %ld(S=%Qd)\n",i,tmpq));
        }
        for (i = 0; i < n; i++)
        {
            OK(SPEX_mpq_mul(tmpq, SPEX_2D(S, 2, i), SPEX_2D(S, 3, i)));

            for (p = 0; p < U->v[i]->nz; p++)
            {
                j = U->v[i]->i[p];
                OK(SPEX_mpz_divexact(tmpz, U->v[i]->x[p],
                    SPEX_MPQ_DEN(tmpq)));
                OK(SPEX_mpz_mul(tmpz, tmpz, SPEX_MPQ_NUM(tmpq)));
                OK(SPEX_gmp_printf("(%ld)%Zd ",j,tmpz));
            }
            OK(SPEX_gmp_printf("......col %ld(S=%Qd)\n",i,tmpq));
        }
        printf("P=[%ld %ld %ld %ld]\n",P[0], P[1], P[2], P[3]);
        printf("Q=[%ld %ld %ld %ld]\n",Q[0], Q[1], Q[2], Q[3]);*/
    }

    return 0;
}

