//------------------------------------------------------------------------------
// SPEX_CHOLMOD/Tcov/tcov_test.c: test coverage for SPEX_CHOLMOD
//------------------------------------------------------------------------------

// SPEX_CHOLMOD: (c) 2020-2021, Jinhao Chen, Timothy A. Davis,
// Erick Moreno-Centeno, Texas A&M University.  All Rights Reserved.  See
// SPEX_CHOLMOD/License for the license.

//------------------------------------------------------------------------------

/*
 * test coverage
 */

#define SPEX_FREE_ALL                            \
{                                                \
    SPEX_FREE(option);                           \
    SPEX_finalize() ;                            \
}

#include "tcov_malloc_test.h"

#define TEST_CHECK(method)                       \
{                                                \
    info = (method) ;                            \
    if (info != SPEX_OK)                         \
    {                                            \
        SPEX_PRINT_INFO (info) ;                 \
        SPEX_FREE_ALL;                           \
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
    //------------------------------------------------------------------
    // Initialize SPEX CHOLMOD process
    //------------------------------------------------------------------

    SPEX_initialize_expert (tcov_malloc, tcov_calloc,
	tcov_realloc, tcov_free) ;

    info = SPEX_initialize ( ) ;
    assert (info == SPEX_PANIC) ;

    //------------------------------------------------------------------
    // Allocate memory
    //------------------------------------------------------------------

    int64_t n=4, i, j, p;
    SPEX_options* option = SPEX_create_default_options();
    if (!option) return 0;//{continue;}

    SPEX_matrix *L = NULL, *U = NULL, *A = NULL;
    SPEX_vector *vk = NULL;
    mpz_t *d = NULL, *sd = NULL;
    mpq_t *S = NULL;
    int64_t *P = NULL, *P_inv = NULL, *Q = NULL, *Q_inv = NULL;

    TEST_CHECK(SPEX_matrix_alloc(&L, 4, 4, true));
    TEST_CHECK(SPEX_matrix_alloc(&U, 4, 4, true));
    int64_t nz = 0;
    for (i = 0; i < 4; i++)
    {
        TEST_CHECK(SPEX_vector_realloc(L->v[i], 4-i));
        TEST_CHECK(SPEX_vector_realloc(U->v[i], 4-i));
        for (p = 0; p < 4-i; p++)
        {
            TEST_CHECK(SPEX_mpz_set_si(L->v[i]->x[p], Lold[nz]));
            TEST_CHECK(SPEX_gmp_printf("L(%ld,%ld)=%Zd=%ld\n", i+p, i,
                L->v[i]->x[p],Lold[nz]));
            TEST_CHECK(SPEX_mpz_set_si(U->v[i]->x[p], Uold[nz]));
            U->v[i]->i[p] = i+p;
            L->v[i]->i[p] = i+p;
            nz++;
        }
        L->v[i]->nz = 4-i;
        U->v[i]->nz = 4-i;
    }
    TEST_CHECK(SPEX_matrix_alloc(&A, 4, 4, true));
    for (i = 0; i < 4; i++)
    {
        TEST_CHECK(SPEX_vector_realloc(A->v[i], 4));
        for (p = 0; p < 4; p++)
        {
            TEST_CHECK(SPEX_mpz_set_si(A->v[i]->x[p], Aold[4*i+p]));
            A->v[i]->i[p] = p;
        }
        A->v[i]->nz = 4;
    }

    TEST_CHECK(SPEX_vector_alloc(&vk, 4, true));
    for (i = 0; i < 4; i++)
    {
        TEST_CHECK(SPEX_mpz_set_ui(vk->x[i], Ak_new[i]));
        vk->i[i] = i;
    }
    vk->nz = 4;

    d  = spex_create_mpz_array(n);
    sd = spex_create_mpz_array(n);
    S  = spex_create_mpq_array(3*n);
    P     = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    Q     = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    P_inv = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    Q_inv = (int64_t*) SPEX_malloc(n*sizeof(int64_t));
    if (!d || !sd || !S || !P || !Q || !P_inv || !Q_inv)
    {
        SPEX_FREE_ALL;
        return 0;
    }
    for (i = 0; i < n; i++)
    {
        TEST_CHECK(SPEX_mpz_set(d[i],  L->v[i]->x[0]));
        TEST_CHECK(SPEX_mpz_set(sd[i], L->v[i]->x[0]));
        TEST_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 1, i), 1, 1));
        TEST_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 2, i), 1, 1));
        TEST_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 3, i), 1, 1));
        P[i] = i;
        Q[i] = i;
        P_inv[i] = i;
        Q_inv[i] = i;
    }
    GOTCHA;
    TEST_CHECK(SPEX_LUU(A, L, U, d, sd, S, P, P_inv, Q, Q_inv, &vk, 1, NULL));

    mpz_t tmpz; SPEX_MPZ_SET_NULL(tmpz);
    TEST_CHECK(SPEX_mpz_init(tmpz));
    mpq_t tmpq; SPEX_MPQ_SET_NULL(tmpq);
    TEST_CHECK(SPEX_mpq_init(tmpq));
    for (i = 0; i < n; i++)
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
    for (i = 0; i < n; i++)
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
    printf("Q=[%ld %ld %ld %ld]\n",Q[0], Q[1], Q[2], Q[3]);


    // reset A=L=U=P=Q=I
    for (i = 0; i < n; i++)
    {
        TEST_CHECK(SPEX_mpz_set_ui(A->v[i]->x[0], 1));
        A->v[i]->i[0] = i;
        A->v[i]->nz = 1;
        TEST_CHECK(SPEX_mpz_set_ui(L->v[i]->x[0], 1));
        L->v[i]->i[0] = i;
        L->v[i]->nz = 1;
        TEST_CHECK(SPEX_mpz_set_ui(U->v[i]->x[0], 1));
        U->v[i]->i[0] = i;
        U->v[i]->nz = 1;
        TEST_CHECK(SPEX_mpz_set(d[i],  L->v[i]->x[0]));
        TEST_CHECK(SPEX_mpz_set(sd[i], L->v[i]->x[0]));
        TEST_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 1, i), 1, 1));
        TEST_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 2, i), 1, 1));
        TEST_CHECK(SPEX_mpq_set_ui(SPEX_2D(S, 3, i), 1, 1));
        P[i] = i;
        Q[i] = i;
        P_inv[i] = i;
        Q_inv[i] = i;
    }
    TEST_CHECK(SPEX_mpz_set_ui(vk->x[0], 2));
    vk->i[0] = 1;
    TEST_CHECK(SPEX_mpz_set_ui(vk->x[1], 5));
    vk->i[1] = 3;
    vk->nz = 2;
    TEST_CHECK(SPEX_LUU(A, L, U, d, sd, S, P, P_inv, Q, Q_inv, &vk, 1, NULL));
    for (i = 0; i < n; i++)
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
    for (i = 0; i < n; i++)
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
    printf("Q=[%ld %ld %ld %ld]\n",Q[0], Q[1], Q[2], Q[3]);

    return 0;
}

