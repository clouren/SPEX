//------------------------------------------------------------------------------
// SPEX_QR/SPEX_QR.h: user #include file for SPEX_QR
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2020, Chris Lourenco, United States Naval Academy. 
// All Rights Reserved.  See SPEX_Cholesky/License for the license.

//------------------------------------------------------------------------------

#ifndef SPEX_QR_H
#define SPEX_QR_H


// This software performs an exact integer-preserving QR factorization
// WARNING: This code is experimental and developmental, please do not use it.


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Authors----------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Christopher Lourenco

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Contact Information----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    Please contact Chris Lourenco (chrisjlourenco@gmail.com, lourenco@usna.edu)
//    or Tim Davis (timdavis@aldenmath.com, DrTimothyAldenDavis@gmail.com,
//                  davis@tamu.edu)


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Copyright--------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    SPEX_QR is free software; you can redistribute it and/or modify
//     it under the terms of either:
//
//        * the GNU Lesser General Public License as published by the
//          Free Software Foundation; either version 3 of the License,
//          or (at your option) any later version.
//
//     or
//
//        * the GNU General Public License as published by the Free Software
//          Foundation; either version 2 of the License, or (at your option) any
//          later version.
//
//    or both in parallel, as here.
//
//    See license.txt for license info.
//
// This software is copyright by Christopher Lourenco
//

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------------DISCLAIMER-----------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// SPEX QR is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//--------------------------Summary---------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

//    To be done


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//---------------------Include files required by SPEX QR------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>
#include "SPEX_Util.h"

// SuiteSparse headers
#include "SuiteSparse_config.h"
#include "colamd.h"
#include "amd.h"

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//-------------------------Default Parameters-----------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// Current version of the code
#define SPEX_QR_VERSION "0.0.1"
#define SPEX_CHOL_VERSION_MAJOR 0
#define SPEX_CHOL_VERSION_MINOR 0
#define SPEX_CHOL_VERSION_SUB   1

#define ASSERT assert

/* Compute the dot product of two integer vectors x,y and return in z */
void SPEX_dot
(
    SPEX_matrix* x,
    SPEX_matrix* y,
    mpz_t z
)
{
    // Check inputs, x and y must both be column vectors of identical size and
    // stored as dense
    ASSERT( x->n == y->n);
    ASSERT( x->m == y->m);
    ASSERT( x->type == SPEX_MPZ);
    ASSERT( y->type == SPEX_MPZ);
    ASSERT( x->kind == SPEX_DENSE);
    ASSERT( y->kind == SPEX_DENSE);
    
    int64_t k;
    
    // Set z = 0
    SPEX_mpz_set_ui(z, 0);
    for (k = 0; k < x->m; k++)
    {
        SPEX_mpz_addmul(z, x->x.mpz[k], y->x.mpz[k]);
    }
}

/* Purpose: Given a matrix A in m*n and B in m*n, compute the dot product of
 * A(:,i) and B(:,j). Assumed to be dense. prod = A(:,i) dot B(:,j)
 */
void SPEX_dense_mat_dot
(
    SPEX_matrix* A,
    int64_t i,
    SPEX_matrix* B,
    int64_t j,
    mpz_t prod
)
{
    ASSERT(A->m == B->m);
    for (int64_t k = 0; k < A->m; k++)
    {
        SPEX_mpz_addmul(prod, SPEX_2D(A, k, i, mpz),
                        SPEX_2D(B, k, j, mpz));
    }
}
    
/* Perform the IPGE version of SPEX QR (aka Algorithm 1 from workpage)
 */
void SPEX_QR_IPGE
(
    SPEX_matrix *A,            // Matrix to be factored
    SPEX_matrix **R_handle,    // upper triangular matrix
    SPEX_matrix **Q_handle     // orthogonal triangular matrix
)
{
    int64_t m = A->m, n = A->n;
    ASSERT( m >= n); // A should be transposed if not true
    ASSERT( A != NULL);
    // Only dense for now
    ASSERT( A->type == SPEX_MPZ);
    ASSERT( A->kind == SPEX_DENSE);
    
    // Indices
    int64_t i, j, k;
    
    // Final matrices Q and R
    SPEX_matrix *Q, *R;
    
    // Allocate R
    SPEX_matrix_allocate(&R, SPEX_DENSE, SPEX_MPZ, n, n, n*n,
        false, true, NULL);
    
    // Set Q = A
     SPEX_matrix_copy(&Q, SPEX_DENSE, SPEX_MPZ, A, NULL);
     
     // Perform Factorization
    for (k = 0; k < n; k++)
    {
        // Compute row k of R
        for (j = k; j < n; j++)
        {
            // R(k,j) = Q(:,k) dot A(:,j)
            // This is very easily parallelized
            SPEX_dense_mat_dot(Q, k, A, j, SPEX_2D(R,k,j,mpz));
        }
        
        // IPGE update Q
        for (i = k+1; i < n; i++)
        {
            for ( j = 0; j < m; j++)
            {
                // Q(j,i) = Q(j,i)*R(k,k)
                SPEX_mpz_mul( SPEX_2D(Q, j, i, mpz),
                              SPEX_2D(Q, j, i, mpz),
                              SPEX_2D(R, k, k, mpz));
                // Q(j,i) = Q(j,i) - R(k,i)*Q(j,k)
                SPEX_mpz_submul( SPEX_2D(Q, j, i, mpz),
                                 SPEX_2D(R, k, i, mpz),
                                 SPEX_2D(Q, j, k, mpz));
                if (k > 0)
                {
                    // Q(j,i) = Q(j,i)/R(k-1,k-1)
                    SPEX_mpz_divexact( SPEX_2D(Q, j, i, mpz),
                                       SPEX_2D(Q, j, i, mpz),
                                       SPEX_2D(R, k-1, k-1, mpz));
                }
            }
        }
    }
    
    (*Q_handle) = Q;
    (*R_handle) = R;
}
                                 
    
    
    
/* Perform the IPGE version of SPEX QR using Pursell method
 */
void SPEX_QR_PURSELL
(
    SPEX_matrix *A,            // Matrix to be factored
    SPEX_matrix **R_handle,    // upper triangular matrix
    SPEX_matrix **Q_handle     // orthogonal triangular matrix
)
{
    
    int64_t m = A->m, n = A->n;
    ASSERT( m >= n); // A should be transposed if not true
    ASSERT( A != NULL);
    // Only dense for now
    ASSERT( A->type == SPEX_MPZ);
    ASSERT( A->kind == SPEX_DENSE);
    
    // Indices
    int64_t i, j, k;
    
    // A_transpose, will be overwritten with Q
    SPEX_matrix *A_T;
    // A2 = A_T A // Will be overwritten with R
    SPEX_matrix *A2;
    
    // Allocate A_T
    SPEX_matrix_allocate(&A_T, SPEX_DENSE, SPEX_MPZ, n, m, n*m,
        false, true, NULL);
    
    // Allocate A2
    SPEX_matrix_allocate(&A2, SPEX_DENSE, SPEX_MPZ, n, n, n*n,
        false, true, NULL);
    
    // Compute A_T
    
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            // A'(i,j) = A(j,i)
            SPEX_mpz_set( SPEX_2D(A_T, i, j, mpz),
                          SPEX_2D(A,   j, i, mpz));
        }
    }
    
    // Compute A2 = A'*A
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            // Compute A2(i,j) = A_T(i,:) dot A(:,j)
            for (k = 0; k < m; k++)
            {
                SPEX_mpz_addmul( SPEX_2D(A2,  i, j, mpz),
                                 SPEX_2D(A_T, i, k, mpz),
                                 SPEX_2D(A,   k, j, mpz));
            }
        }
    }
    
    // Now, the factorization is computed by performing IPGE on [A2 AT]
    // We can either construct this matrix directly or do it more efficiently
    
    // Compute R and Q directly
    // A2 = A'*A will be overwriten with R (dimension n*n)
    // A' will be overwritten with Q (dimension n*m)
    
    // Perform IPGE on A2
    for (k = 0; k < n; k++)
    {
        for (i = k+1; i < n; i++)
        {
            for (j = k+1; j < n; j++)
            {
                if (i > j)
                {
                    // i > j -> entries in lower triangular portion, set them to 0
                    SPEX_mpz_set_ui( SPEX_2D(A2, i, j, mpz), 0);
                }
                else
                {
                    // A2(i,j) = A2(i,j)*A2(k,k)
                    SPEX_mpz_mul ( SPEX_2D(A2, i, j, mpz),
                                   SPEX_2D(A2, i, j, mpz),
                                   SPEX_2D(A2, k, k, mpz));
                    // A2(i,j) = A2(i,j) - A2(i,k)*A2(i,k)
                    SPEX_mpz_submul( SPEX_2D(A2, i, j, mpz),
                                     SPEX_2D(A2, k, i, mpz),
                                     SPEX_2D(A2, k, j, mpz));
                    if (k > 0)
                    {
                        SPEX_mpz_divexact( SPEX_2D(A2, i, j, mpz),
                                           SPEX_2D(A2, i, j, mpz),
                                           SPEX_2D(A2, k-1, k-1, mpz));
                    }
                }
            }
        }
        // Now that we have A2, we can compute Q
        for (i = k+1; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                SPEX_mpz_mul( SPEX_2D(A_T, i, j, mpz),
                              SPEX_2D(A2,  k, k, mpz),
                              SPEX_2D(A_T, i, j, mpz));
                SPEX_mpz_submul( SPEX_2D(A_T, i, j, mpz),
                                 SPEX_2D(A2,  k, i, mpz),
                                 SPEX_2D(A_T, k, j, mpz));
                if (k > 0)
                {
                    SPEX_mpz_divexact( SPEX_2D(A_T, i, j, mpz),
                                       SPEX_2D(A_T, i, j, mpz),
                                       SPEX_2D(A2, k-1, k-1, mpz));
                }
            }
        }
    }
    
    for (k = 1; k < n; k++)
        SPEX_mpz_set_ui( SPEX_2D(A2, k, 0, mpz), 0);
                                    
    (*Q_handle) = A_T;
    (*R_handle) = A2;
    
}
 
 
/* Perform the IPGE version of SPEX QR using Pursell method
 */
void SPEX_QR_PURSELL2
(
    SPEX_matrix *A,            // Matrix to be factored
    SPEX_matrix **R_handle,    // upper triangular matrix
    SPEX_matrix **Q_handle     // orthogonal triangular matrix
)
{
    
    int64_t m = A->m, n = A->n;
    ASSERT( m >= n); // A should be transposed if not true
    ASSERT( A != NULL);
    // Only dense for now
    ASSERT( A->type == SPEX_MPZ);
    ASSERT( A->kind == SPEX_DENSE);
    
    // Indices
    int64_t i, j, k;
    
    // Final matrices Q and R
    SPEX_matrix *Q, *R;
    
    // A_transpose
    SPEX_matrix *A_T;
    // A2 = A_T A
    SPEX_matrix *A2;
    
    // Allocate R
    SPEX_matrix_allocate(&R, SPEX_DENSE, SPEX_MPZ, n, n, n*n,
        false, true, NULL);
    
    // Allocate Q
    SPEX_matrix_allocate(&Q, SPEX_DENSE, SPEX_MPZ, m, n, m*n,
        false, true, NULL);
     
    // Allocate A_T
    SPEX_matrix_allocate(&A_T, SPEX_DENSE, SPEX_MPZ, n, m, n*m,
        false, true, NULL);
    
    // Allocate A2
    SPEX_matrix_allocate(&A2, SPEX_DENSE, SPEX_MPZ, n, n, n*n,
        false, true, NULL);
    
    // Compute A_T
    
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            SPEX_mpz_set( SPEX_2D(A_T, i, j, mpz),
                          SPEX_2D(A, j, i, mpz));
        }
    }
    
    // Compute A2 = A_T*A
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            // Compute A2(i,j) = A_T(i,:) dot A(:,j)
            for (k = 0; k < m; k++)
            {
                SPEX_mpz_addmul( SPEX_2D(A2, i, j, mpz),
                                 SPEX_2D(A_T, i, k, mpz),
                                 SPEX_2D(A, k, j, mpz));
            }
        }
    }
    
    // Now, the factorization is computed by performing IPGE on [A2 AT]
    // We can either construct this matrix directly or do it more efficiently
    
    // This is the A3 version
    SPEX_matrix *A3;
    // Allocate A3
    SPEX_matrix_allocate(&A3, SPEX_DENSE, SPEX_MPZ, n, n+m, n*(n+m),
        false, true, NULL);
    // Populate first n*n portion
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            SPEX_mpz_set( SPEX_2D(A3, i, j, mpz),
                          SPEX_2D(A2, i, j, mpz));
        }
    }
    
    // Populate last n*m portion
    for(i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            SPEX_mpz_set( SPEX_2D(A3, i, j+n, mpz),
                          SPEX_2D(A_T, i, j, mpz));
        }
    }
    
    // Now, perform IPGE on A3
    for (k = 0; k < n; k++)
    {
        for (i = k+1; i < n; i++)
        {
            for (j = k+1; j < A3->n; j++)
            {
                SPEX_mpz_mul ( SPEX_2D(A3, i, j, mpz),
                               SPEX_2D(A3, i, j, mpz),
                               SPEX_2D(A3, k, k, mpz));
                SPEX_mpz_submul( SPEX_2D(A3, i, j, mpz),
                                 SPEX_2D(A3, i, k, mpz),
                                 SPEX_2D(A3, k, j, mpz));
                if (k > 0)
                    SPEX_mpz_divexact( SPEX_2D(A3, i, j, mpz),
                                       SPEX_2D(A3, i, j, mpz),
                                       SPEX_2D(A3, k-1, k-1, mpz));
            }
        }
    }
    
    // Now get Q and R
    // Populate first n*n portion
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (i <= j)
            {
            SPEX_mpz_set( SPEX_2D(R, i, j, mpz),
                          SPEX_2D(A3, i, j, mpz));
            }
            else
                SPEX_mpz_set_ui( SPEX_2D(R,i,j,mpz), 0);
        }
    }
    
    // Populate last n*m portion
    for(i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            SPEX_mpz_set( SPEX_2D(Q, i, j, mpz),
                          SPEX_2D(A3, i, j+n, mpz));
        }
    }
    
    (*Q_handle) = Q;
    (*R_handle) = R;
}
    
    

#endif
