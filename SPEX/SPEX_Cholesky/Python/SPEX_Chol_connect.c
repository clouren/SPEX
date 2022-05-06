//------------------------------------------------------------------------------
// SPEX_Cholesky/Python/SPEX_Chol_connect.c: use SPEX_Chol in Python
//------------------------------------------------------------------------------

// SPEX_Cholesky: (c) 2021, Chris Lourenco, United States Naval Academy, 
// Lorena Mejia Domenzain, Erick Moreno-Centeno, Timothy A. Davis,
// Texas A&M University. All Rights Reserved. 
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

#include "SPEX_Chol_connect.h"

#define FREE_WORKSPACE              \
    SPEX_matrix_free(&A, option);   \
    SPEX_matrix_free(&b, option);   \
    SPEX_matrix_free(&x, option);   \
    SPEX_matrix_free(&A_in, option); \
    SPEX_matrix_free(&b_in, option);   \
    SPEX_FREE(option);              \
    SPEX_finalize();                

SPEX_info SPEX_python_backslash
( 
//output
 char** sol_char,      // solution (null if double output)
 double* sol_doub,     // solution (null if char output)
 void** sol_void,
 //input
 bool charOut,    // True if char** output, false if double
 int64_t* Ap,     // column pointers of A, an array size is n+1 
 int64_t* Ai,     // row indices of A, of size nzmax.
 double* Ax,      // values of A
 double* bx,      // values of b
 int m,           // Number of rows of A
 int n,           // Number of columns of A
 int nz,           // Number of nonzeros in A
 int ordering     // type of ordering: 0-none, 1-colamd, 2-amd
)
{
    //this function will make the SPEX_matrix A, b and x and call SPEX_Chol_backslash
    SPEX_info info;

    //--------------------------------------------------------------------------
    // Prior to using SPEX Chol, its environment must be initialized. This is
    // done by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------
    SPEX_initialize();
    
    //--------------------------------------------------------------------------
    // Validate Inputs
    //--------------------------------------------------------------------------

    if (!Ap || !Ai || !Ax || !bx)
    {
        return SPEX_INCORRECT_INPUT;
    }

    // ordering must be acceptable
    if (ordering != 0 && ordering != 1 && ordering != 2)
    {
        return SPEX_INCORRECT_INPUT;
    }

    if (n == 0 || m == 0 || n != m)
    {
        return SPEX_INCORRECT_INPUT;
    }

    //--------------------------------------------------------------------------
    // Declare our data structures
    //--------------------------------------------------------------------------
    SPEX_matrix* A_in = NULL;       //input matrix
    SPEX_matrix* b_in = NULL;       //input rhs
    SPEX_matrix* A = NULL;          //copy of input matrix in CSC MPZ 
    SPEX_matrix* b = NULL;          //copy of input rhs in CSC MPZ 
    SPEX_matrix* x = NULL;          //solution
    
    SPEX_options *option = NULL;
    SPEX_create_default_options(&option); //TOFREE opitons
    SPEX_preorder order_in = ordering;
    option->order = ordering; 

    //--------------------------------------------------------------------------
    // Allocate memory, populate in A and b
    //--------------------------------------------------------------------------
    SPEX_matrix_allocate(&A_in, SPEX_CSC, SPEX_FP64, n, n, nz, true, true, option); 
    SPEX_matrix_allocate(&b_in, SPEX_DENSE, SPEX_FP64, n, 1, n, true, true, option);

    A_in->p=Ap;
    A_in->i=Ai;
    A_in->x.fp64=Ax;

    b_in->x.fp64=bx;

    // At this point, A_in is a shallow double CSC matrix. We make a copy of it with A
    // A is a CSC matrix with mpz entries
    SPEX_matrix_copy(&A, SPEX_CSC, SPEX_MPZ, A_in, option);
    // b is a dense matrix with mpz entries
    SPEX_matrix_copy(&b, SPEX_DENSE, SPEX_MPZ, b_in, option);

    //--------------------------------------------------------------------------
    // solve
    //-------------------------------------------------------------------------
    //option->check = true;
    //option->print_level = 2;
    
    SPEX_Chol_backslash(&x, SPEX_MPQ, A, b, option); 

    //--------------------------------------------------------------------------
    // Return output as desired type
    //-----------------------------------------------------------------------
    if(charOut)
    {
        for (int i = 0; i < n; ++i)
        {
            char *s ;
            int status = SPEX_mpfr_asprintf (&s, "%Qd", x->x.mpq [i]) ;
            if (status < 0)
            {
                printf("error converting x to string") ;
            }
            //check string size
            int sizeStr;
            sizeStr=strlen(s);
            //allocate sol_char[i]
            sol_char[i] = (char*) malloc (sizeStr*sizeof(char));
            //copy s into sol_char[i]
            strcpy(sol_char[i],s);
        }
        for (int i = 0; i < n; ++i)
        {
            char *s ;
            int status = SPEX_mpfr_asprintf (&s, "%Qd", x->x.mpq [i]) ;
            if (status < 0)
            {
                printf("error converting x to string") ;
            }
            //check string size
            int sizeStr;
            sizeStr=strlen(s);
            //allocate sol_char[i]
            sol_void[i] = (void*) malloc (sizeStr*sizeof(char));
            //copy s into sol_char[i]
            strcpy(sol_void[i],s);
        }
    }
    else
    {
        SPEX_matrix* x2 = NULL ;
        SPEX_matrix_copy(&x2, SPEX_DENSE, SPEX_FP64, x, option);
        for (int i = 0; i < n; ++i)
        {   
            sol_doub[i]=x2->x.fp64[i];
        }
        for (int i = 0; i < n; ++i)
        {
            //sol_void[i] = (void*) malloc (1*sizeof(double));
            (sol_void[i])=&(x2->x.fp64[i]);
        }  
    }
    printf("\n void mem %p,%p,%p\n",sol_void[0],sol_void[1],sol_void[2]);
    
    FREE_WORKSPACE;
}



void python_backslash_double
( 
 //output
 double* sol,                // solution
 //input
 int64_t* col_pointers,     // column pointers of A, an array size is n+1
 int64_t* row_index,        // row indices of A, of size nzmax.
 double* data,              // values of A
 double* rhs,               // values of b
 int n,                     // # of rows/columns of A
 int nz,                    // max # entries for CSC matrix A
 int ordering               // type of ordering: 0-none, 1-colamd, 2-amd
)
{
    //this function will make the SPEX_matrix A, b and x and call SPEX_Chol_backslash

    //--------------------------------------------------------------------------
    // Prior to using SPEX Chol, its environment must be initialized. This is
    // done by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------
    SPEX_initialize();
    /**/
    //--------------------------------------------------------------------------
    // Declare our data structures
    //--------------------------------------------------------------------------
    SPEX_matrix* A_in = NULL;       //input matrix
    SPEX_matrix* b_in = NULL;       //input rhs
    SPEX_matrix* A = NULL;          //copy of input matrix in CSC MPZ 
    SPEX_matrix* b = NULL;          //copy of input rhs in CSC MPZ 
    SPEX_matrix* x = NULL;          //solution
    
    SPEX_options *option = NULL;
    SPEX_create_default_options(&option); 
    SPEX_preorder order_in = ordering;
    option->order = ordering; 
    
    //--------------------------------------------------------------------------
    // Allocate memory, populate in A and b
    //--------------------------------------------------------------------------
    SPEX_matrix_allocate(&A_in, SPEX_CSC, SPEX_FP64, n, n, nz, false, true, option); 
    SPEX_matrix_allocate(&b_in, SPEX_DENSE, SPEX_FP64, n, 1, n, false, true, option);

    
    //populate A_in and b_in with function inputs 
    for (int i = 0; i < n; ++i)
    {
        A_in->p[i]=col_pointers[i];
        b_in->x.fp64[i] = rhs[i];
    }
    A_in->p[n]=col_pointers[n];
    for (int i = 0; i < nz; ++i)
    {
        A_in->i[i]=row_index[i];
        A_in->x.fp64[i]=data[i];
    }
    
    // At this point, A_in is a double CSC matrix. We make a copy of it with A
    // A is a CSC matrix with mpz entries
    SPEX_matrix_copy(&A, SPEX_CSC, SPEX_MPZ, A_in, option);
    // b is a dense matrix with mpz entries
    SPEX_matrix_copy(&b, SPEX_DENSE, SPEX_MPZ, b_in, option);

    //--------------------------------------------------------------------------
    // solve
    //-------------------------------------------------------------------------
    //option->check = true;
    //option->print_level = 2;
    
    SPEX_Chol_backslash(&x, SPEX_FP64, A, b, option); 
    
    
    for (int i = 0; i < n; ++i)
    {
        sol[i]=x->x.fp64[i];
    }

    FREE_WORKSPACE;/**/
}

void python_backslash_char
( 
 //output
 char** sol,                //solution
 //input
 int64_t* col_pointers,     //column pointers of A, an array size is n+1
 int64_t* row_index,        //row indices of A, of size nzmax.
 double* data,              //values of A
 double* rhs,               //values of b
 int n,                     // # of rows/columns of A
 int nz,                    // max # entries for CSC matrix A
 int ordering               // type of ordering: 0-none, 1-colamd, 2-amd
)
{
    //basically the same funcion as the previous one but with different solution type
    //this function will make the SPEX_matrix A, b and x and call SPEX_Chol_backslash

    //--------------------------------------------------------------------------
    // Prior to using SPEX Chol, its environment must be initialized. This is
    // done by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------
    SPEX_initialize();
    /**/
    //--------------------------------------------------------------------------
    // Declare our data structures
    //--------------------------------------------------------------------------
    SPEX_matrix* A_in = NULL;       //input matrix
    SPEX_matrix* b_in = NULL;       //input rhs
    SPEX_matrix* A = NULL;          //copy of input matrix in CSC MPZ 
    SPEX_matrix* b = NULL;          //copy of input rhs in CSC MPZ 
    SPEX_matrix* x = NULL;          //solution
    
    SPEX_options *option = NULL;
    SPEX_create_default_options(&option); //TOFREE opitons
    
    //--------------------------------------------------------------------------
    // Allocate memory, populate in A and b
    //--------------------------------------------------------------------------
    SPEX_matrix_allocate(&A_in, SPEX_CSC, SPEX_FP64, n, n, nz, false, true, option); 
    SPEX_matrix_allocate(&b_in, SPEX_DENSE, SPEX_FP64, n, 1, n, false, true, option);

    
    //populate A_in and b_in with function inputs 
    for (int i = 0; i < n; ++i)
    {
        A_in->p[i]=col_pointers[i];
        b_in->x.fp64[i] = rhs[i];
    }
    A_in->p[n]=col_pointers[n];
    for (int i = 0; i < nz; ++i)
    {
        A_in->i[i]=row_index[i];
        A_in->x.fp64[i]=data[i];
    }

    // At this point, A_in is a double CSC matrix. We make a copy of it with A
    // A is a CSC matrix with mpz entries
    SPEX_matrix_copy(&A, SPEX_CSC, SPEX_MPZ, A_in, option);
    // b is a dense matrix with mpz entries
    SPEX_matrix_copy(&b, SPEX_DENSE, SPEX_MPZ, b_in, option);

    //--------------------------------------------------------------------------
    // solve
    //-------------------------------------------------------------------------
    //option->check = true;
    //option->print_level = 2;
    
    SPEX_Chol_backslash(&x, SPEX_MPQ, A, b, option); 

    //allocate sol
    char** sol2;
    sol2 = (char**) malloc (n*sizeof(char*));
    for (int i = 0; i < n; ++i)
    {
        //hacerlo string
        char *s ;
        int status = SPEX_mpfr_asprintf (&s, "%Qd", x->x.mpq [i]) ;
        if (status < 0)
        {
            printf("error converting x to string") ;
        }
        //checar tamanio
        int sizeStr;
        sizeStr=strlen(s);
        //falta hacer allocation de sol[i]
        sol2[i] = (char*) malloc (sizeStr*sizeof(char));
        //copiarlo en sol[i]
        strcpy(sol2[i],s);
        sol[i]=sol2[i];

    }
    //printf("%s,%s,%s\n",sol2[0],sol2[1],sol2[2]);
    //posible solucion, regresar todo como una gran cadena numer/denom,numer/denom,numer/denom
    
    FREE_WORKSPACE;/**/
}


