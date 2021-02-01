//------------------------------------------------------------------------------
// SPEX_WAMF/SPEX_WAMF_get_weights: Obtain the weights of a matrix
//------------------------------------------------------------------------------

// SPEX_WAMF: (c) 2021, Chris Lourenco US Naval Academy, Erick Moreno-Centeno, Texas
// A&M University. All Rights Reserved.  See SPEX_WAMF/License for the license.

//------------------------------------------------------------------------------

#include "SPEX_WAMF.h"

/* Purpose: WAMF stores a weight associated with each potential pivot row/column
 * (similar to the degree in minimum degree). This function initializes the weights
 * by creating an initial set of weights of the graph. 
 * 
 * Input arguments:
 * 
 *  A: matrix to be analyzed
 * 
 *  option: an int that determines how the weights are calculated. Currently these are:
 *      0: w[k] = sum A(i,k)*i
 *      1: w[k] = sum A(:,k)
 *      2: w[k] = mean (A(:,k))
 *      3: w[k] = min A(:,k)
 *      4: w[k] = max A(:,k)
 *      5: w[k] = A(k,k), if A(k,k) = 0 sum.
 * 
 * 
 * On success, this function returns a double* array containing the weights. On failure 
 * it returns NULL
 */

double* SPEX_WAMF_get_weights
(
    SPEX_matrix* A,
    int64_t option
)
{
    if ( !A || !A->x.fp64 || !A->i || !A->p || option < 0 || option > 5) {return NULL;}
    ASSERT(A->kind == SPEX_CSC);
    ASSERT(A->type == SPEX_FP64);
    
    int64_t k, i, p, nnz, n = A->n;
    double* w = SPEX_calloc(n, sizeof(double));
    
    if (w == NULL) return NULL;
    
       
    // Weighted sum: w[k] =  sum (A(i,k)*i)
    if (option == 0)
    {
        for (k = 0; k < n; k++)
        {
            for (i = A->p[k]; i < A->p[k+1]; i++)
            {
                w[k] += A->x.fp64[i] * A->i[i];
            }
        }
    }
    
    // w[k] = sum(A(:,k))
    else if (option == 1)
    {
        for (k = 0; k < n; k++)
        {
            for (i = A->p[k]; i < A->p[k+1]; i++)
            {
                w[k] += A->x.fp64[i];
            }
        }
    }
    
    // w[k] = mean(A(:,k))
    else if (option == 2)
    {
        for (k = 0; k < n; k++)
        {
            for (i = A->p[k]; i < A->p[k+1]; i++)
            {
                w[k] += A->x.fp64[i];
            }
            nnz = A->p[k+1] - A->p[k];
            w[k] = (double) w[k] / nnz;
            w[k] = ceil(w[k]);
        }
    }
    
    // w[k] = min(A(:,k))
    else if (option == 3)
    {
        for (k = 0; k < n; k++)
        {
            w[k] = A->x.fp64[A->p[k]];
            for (i = A->p[k]+1; i < A->p[k+1]; i++)
            {
                if (w[k] > A->x.fp64[i])
                {
                    w[k] = A->x.fp64[i];
                }
            }
        }
    }
    
    // w[k] = max(A(:,k))
    else if (option == 4)
    {
        for (k = 0; k < n; k++)
        {
            w[k] = A->x.fp64[A->p[k]];
            for (i = A->p[k]+1; i < A->p[k+1]; i++)
            {
                if (w[k] < A->x.fp64[i])
                {
                    w[k] = A->x.fp64[i];
                }
            }
        }
    }
    
    // w[k] = A(k,k), if A(k,k) = 0, sum
    else if (option == 5)       
    {
        int64_t count = 0;
        for (k = 0; k < n; k++)
        {
            for (i = A->p[k]; i < A->p[k+1]; i++)
            {
                if ( A->i[i] == k)
                {
                    count+=1;
                    w[k] = A->x.fp64[i];
                    break;
                }
                else
                    w[k] += A->x.fp64[i];
            }
        }
    }
    return w;
}
