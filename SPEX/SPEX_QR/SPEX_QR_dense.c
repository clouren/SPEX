
# include "SPEX_QR.h"

int main( int argc, char* argv[] )
{

    //--------------------------------------------------------------------------
    // Prior to using SPEX QR, its environment must be initialized. This is done
    // by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------

    SPEX_initialize();

    //--------------------------------------------------------------------------
    // Declare and initialize essential variables
    //--------------------------------------------------------------------------

    SPEX_info ok;
    int64_t n = 300, nz = n*n, num=0;
    SPEX_matrix *A = NULL ;                     // input matrix
    SPEX_matrix *A2 = NULL ;                    // input matrix (to be generated)
    SPEX_matrix *R = NULL;                      // Upper triangular matrix
    SPEX_matrix *Q = NULL;                      // Orthogonal Matrix
    SPEX_options *option = NULL;
    SPEX_create_default_options(&option);
    if (!option)
    {
        fprintf (stderr, "Error! OUT of MEMORY!\n");
        //FREE_WORKSPACE;
        return 0;
    }

    //--------------------------------------------------------------------------
    // Generate a random dense matrix
    //--------------------------------------------------------------------------

    // A2 is a n*n triplet matrix whose entries are FP64 Note that the first
    // boolean parameter says that the matrix is not shallow, so that A2->i,
    // A2->j, and A2->x are calloc'd. The second boolean parameter is meaningless
    // for FP64 matrices, but it tells SPEX LU to allocate the values of A2->x
    // for the mpz_t, mpq_t, and mpfr_t entries
    SPEX_matrix_allocate(&A2, SPEX_TRIPLET, SPEX_FP64, n, n, nz,
        false, true, option);
    
    // Randomly generate the input
    unsigned int seed = 10;
    srand(seed);
    for (int64_t k = 0; k < n; k++)
    {
        for (int64_t p = 0; p < n; p++)
        {
            A2->i[num] = k;
            A2->j[num] = p;
            //A2->x.fp64[num] = rand(); rdno = rand() % (b-a) + a
            A2->x.fp64[num] = rand() % (10-1) + 1;
            num+=1;
        }
    }

    A2->nz = n*n;
    
    // Create A as a copy of A2
    // A is a copy of the A2 matrix. A is a dense matrix with mpz_t entries
    SPEX_matrix_copy(&A, SPEX_DENSE, SPEX_MPZ, A2, option);
    
    //--------------------------------------------------------------------------
    // Factorize
    //--------------------------------------------------------------------------
    
    // Perform the mulitplication heavy QR_IPGE. This method is focused on doing
    // dot products and tries to limit the number of divisions
    // Better for parallelization and memory
    
    clock_t start_solve1 = clock();
    
    SPEX_QR_IPGE( A, &R, &Q);
    
    clock_t end_solve1 = clock();
    
    // Use IPGE Pursell method. This approach performs IPGE on [A' * A | A']
    // This specific version operates on A'*A and A' seperately and is quite inefficient
    
    SPEX_matrix *Q2, *R2;
    
    clock_t start_solve2 = clock();
    
    SPEX_QR_PURSELL( A, &R2, &Q2);
    
    clock_t end_solve2 = clock();
    
    // Use IPGE Pursell method. This approach explicitly constructs [A'*A | A']
    // The pursell method is generally quite division heavy.
    
    SPEX_matrix *Q3, *R3;
    
    clock_t start_solve3 = clock();
    
    SPEX_QR_PURSELL2( A, &R3, &Q3);
    
    clock_t end_solve3 = clock();
/*    
    printf("\nA is: \n");
    for (int64_t i = 0; i < A->m; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            gmp_printf(" %Zd", SPEX_2D(A,i,j,mpz));
        }
        printf("\n");
    }
    
    printf("\nOur Q is: \n");
    for (int64_t i = 0; i < A->m; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            gmp_printf(" %Zd", SPEX_2D(Q,i,j,mpz));
        }
        printf("\n");
    }
    
    printf("\nOur R is: \n");
    for (int64_t i = 0; i < A->n; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            gmp_printf(" %Zd", SPEX_2D(R,i,j,mpz));
        }
        printf("\n");
    }
    
    printf("\nPursell 1 Q is: \n");
    for (int64_t i = 0; i < A->m; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            gmp_printf(" %Zd", SPEX_2D(Q2,i,j,mpz));
        }
        printf("\n");
    }
    
    printf("\nPursell 1 R is: \n");
    for (int64_t i = 0; i < A->n; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            gmp_printf(" %Zd", SPEX_2D(R2,i,j,mpz));
        }
        printf("\n");
    }
    
    printf("\nPursell 2 Q is: \n");
    for (int64_t i = 0; i < A->m; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            gmp_printf(" %Zd", SPEX_2D(Q3,i,j,mpz));
        }
        printf("\n");
    }
    
    printf("\nPursell 2 R is: \n");
    for (int64_t i = 0; i < A->n; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            gmp_printf(" %Zd", SPEX_2D(R3,i,j,mpz));
        }
        printf("\n");
    }
*/    
    // Now check to make sure they are the same
    for (int64_t i = 0; i < A->m; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            int r ;
            // Have to transpose Q here because theirs is backwards
            SPEX_mpz_cmp(&r, SPEX_2D(Q,j,i,mpz), SPEX_2D(Q2, i, j, mpz));
            if ( r != 0)
                printf("\n Q2 Incorrect at %ld %ld", i, j);
        }
    }
    
    for (int64_t i = 0; i < A->n; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            int r ;
            SPEX_mpz_cmp(&r, SPEX_2D(R,i,j,mpz), SPEX_2D(R2, i, j, mpz));
            if ( r != 0)
                printf("\n R2 Incorrect at %ld %ld", i, j);
        }
    }
    
    for (int64_t i = 0; i < A->m; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            int r ;
            SPEX_mpz_cmp(&r, SPEX_2D(Q,j,i,mpz), SPEX_2D(Q3, i, j, mpz));
            if ( r != 0)
                printf("\n Q3 Incorrect at %ld %ld", i, j);
        }
    }
    
    for (int64_t i = 0; i < A->n; i++)
    {
        for (int64_t j = 0; j < A->n; j++)
        {
            int r ;
            SPEX_mpz_cmp(&r, SPEX_2D(R,i,j,mpz), SPEX_2D(R3, i, j, mpz));
            if ( r != 0)
                printf("\n R3 Incorrect at %ld %ld", i, j);
        }
    }
    
    //--------------------------------------------------------------------------
    // Output & Timing Stats
    //--------------------------------------------------------------------------
    
    double t_qr1 =  (double) (end_solve1 - start_solve1) / CLOCKS_PER_SEC;
    double t_qr2 =  (double) (end_solve2 - start_solve2) / CLOCKS_PER_SEC;
    double t_qr3 =  (double) (end_solve3 - start_solve3) / CLOCKS_PER_SEC;


    double rat  = t_qr2 / t_qr1;
    double rat2 = t_qr3 / t_qr1;
    
    printf("\nIPGE QR time: \t\t\t%lf\n\n", t_qr1);
    printf("\nIPGE QR Pursell time: \t\t%lf\n\n", t_qr2);
    printf("\nIPGE QR Pursell2 time: \t\t%lf\n\n", t_qr3);
    printf("\nRatio IPGE QR/ Pursell1: \t%lf\n\n", rat);
    printf("\nRatio IPGE QR/ Pursell2: \t%lf\n\n", rat2);
    
    //--------------------------------------------------------------------------
    // Free Memory
    //--------------------------------------------------------------------------
    
                 
}
    
