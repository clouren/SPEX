// example.c
#include "example.h"

#define FREE_WORKSPACE              \
    SPEX_matrix_free(&A, option);   \
    SPEX_matrix_free(&b, option);   \
    SPEX_matrix_free(&x, option);   \
    SPEX_matrix_free(&A_in, option); \
    SPEX_matrix_free(&b_in, option);   \
    SPEX_FREE(option);              \
    SPEX_finalize();                

/*
when i use:     
    SPEX_matrix_free(&A_in, option);   \
    SPEX_matrix_free(&b_in, option);   \
 I get "double free or corruption (out)" TODO FIX
*/
//TODO clean file, remove everything but python_backslash2 (change name)
//TODO valgrind for memory related things :)

void read(double* input, int size)
{
  int i;
  for(i=0;i<size;i++)
    input[i] = i;
}


int python_backslash
( 
 //output
 double** sol,                   //solution
 //input
 int64_t* col_pointers,     //column pointers of A, an array size is n+1
 int64_t* row_index,        //row indices of A, of size nzmax.
 double* data,              //values of A
 double* rhs,               //values of b
 int n,                     // # of rows/columns of A
 int nz                     // max # entries for CSC matrix A
)
{
    //eventually make into void function (or something) right now return is to see if things are working
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
    
    //--------------------------------------------------------------------------
    // Allocate memory, populate in A and b
    //--------------------------------------------------------------------------
    SPEX_matrix_allocate(&A_in, SPEX_CSC, SPEX_FP64, n, n, nz, false, true, option);
    SPEX_matrix_allocate(&b_in, SPEX_DENSE, SPEX_FP64, n, 1, n, false, true, option);

    
    //populate A_in with function inputs
    A_in->p = col_pointers;
    (A_in->i) = row_index;
    A_in->x.fp64 = data;
    //Populate b_in with function inputs
    b_in->x.fp64 = rhs;
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
    
    *sol=x->x.fp64;
    //printf("%f,%f,%f\n",x->x.fp64[0],x->x.fp64[1],x->x.fp64[2]);
    //printf("%f,%f,%f\n",sol[0],sol[1],sol[2]);
    // x2 is a copy of the solution. x2 is a dense matrix with double entries (can I do strings???)
    //SPEX_matrix_copy(&x2, SPEX_DENSE, SPEX_FP64, x, option);
    
    //FREE_WORKSPACE;/**/
    return 1;
}

void python_backslash2
( 
 //output
 double* sol,                //solution
 //input
 int64_t* col_pointers,     //column pointers of A, an array size is n+1
 int64_t* row_index,        //row indices of A, of size nzmax.
 double* data,              //values of A
 double* rhs,               //values of b
 int n,                     // # of rows/columns of A
 int nz                     // max # entries for CSC matrix A
)
{
    //eventually make into void function (or something) right now return is to see if things are working
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
    SPEX_matrix_allocate(&A_in, SPEX_CSC, SPEX_FP64, n, n, nz, false, true, option); //TOFREE A_in (x4?)
    SPEX_matrix_allocate(&b_in, SPEX_DENSE, SPEX_FP64, n, 1, n, false, true, option); //TOFREE b_in

    
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
    
    //A_in->p = col_pointers;
    //(A_in->i) = row_index;
    //A_in->x.fp64 = data;
    //Populate b_in with function inputs
    //b_in->x.fp64 = rhs;
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
    //printf("%f,%f,%f\n",x->x.fp64[0],x->x.fp64[1],x->x.fp64[2]);
    //printf("%f,%f,%f\n",sol[0],sol[1],sol[2]);
    // x2 is a copy of the solution. x2 is a dense matrix with double entries (can I do strings???)
    //SPEX_matrix_copy(&x2, SPEX_DENSE, SPEX_FP64, x, option);
    
    FREE_WORKSPACE;/**/
}

