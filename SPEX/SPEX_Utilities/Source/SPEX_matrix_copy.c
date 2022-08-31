//------------------------------------------------------------------------------
// SPEX_Utilities/SPEX_matrix_copy: create a copy of a matrix
//------------------------------------------------------------------------------

// SPEX_Utilities: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------

// SPEX_matrix_copy creates a SPEX_matrix C that is a modified copy of a
// SPEX_matrix A.  The new matrix C can have a different kind and type
// than A.

// The input matrix A is assumed to be valid.  It can be checked first with
// SPEX_matrix_check, if desired.  If the input matrix A is not valid, results
// are undefined.

// SPEX supports 16 matrix formats:  15 of them are all combinations of
// (CSC, triplet, dense) x (mpz, mpq, mpfr, int64, double).  The 16th format
// is dynamic CSC, which can only be mpz.  This function can convert an input
// matrix A in any of these 16 formats, into an output matrix C in any of the
// 16 supported formats.

#define SPEX_FREE_WORK                  \
    SPEX_matrix_free (&T, option) ;     \
    SPEX_matrix_free (&Y, option) ;     \
    SPEX_FREE (W) ;

#define SPEX_FREE_ALL                   \
    SPEX_FREE_WORK ;                    \
    SPEX_matrix_free (&C, option) ;

#include "spex_util_internal.h"

SPEX_info SPEX_matrix_copy
(
    SPEX_matrix **C_handle, // matrix to create (never shallow)
    // inputs, not modified:
    SPEX_kind C_kind,       // C->kind: CSC, triplet, dense, or dynamic_CSC
    SPEX_type C_type,       // C->type: mpz_t, mpq_t, mpfr_t, int64_t, or double
    const SPEX_matrix *A,         // matrix to make a copy of (may be shallow)
    const SPEX_options *option
)
{

    //--------------------------------------------------------------------------
    // check inputs
    //--------------------------------------------------------------------------
    SPEX_info info ;
    if (!spex_initialized ( )) return (SPEX_PANIC) ;

    int64_t nz;
    SPEX_matrix *C = NULL ;
    SPEX_matrix *Y = NULL ;
    SPEX_matrix *T = NULL ;
    int64_t *W = NULL ;

    SPEX_CHECK (SPEX_matrix_nnz (&nz, A, option)) ;
    ASSERT( nz >= 0);
    if (C_handle == NULL || nz < 0 ||
      //checked in SPEX_matrix_nnz
      //A == NULL || A->kind < SPEX_CSC || A->kind > SPEX_DYNAMIC_CSC ||
        A->type < SPEX_MPZ || A->type > SPEX_FP64  ||
        C_kind  < SPEX_CSC || C_kind  > SPEX_DYNAMIC_CSC ||
        C_type  < SPEX_MPZ || C_type  > SPEX_FP64 ||
        (A->kind == SPEX_DYNAMIC_CSC && A->type != SPEX_MPZ) ||
        (C_kind  == SPEX_DYNAMIC_CSC && C_type  != SPEX_MPZ))
    {
        return (SPEX_INCORRECT_INPUT) ;
    }
    (*C_handle) = NULL ;
    int64_t m = A->m ;
    int64_t n = A->n ;
    mpfr_rnd_t round = SPEX_OPTION_ROUND (option) ;

    //--------------------------------------------------------------------------
    // copy and convert A into C
    //--------------------------------------------------------------------------

    switch (C_kind)
    {

        //----------------------------------------------------------------------
        // C is CSC
        //----------------------------------------------------------------------

        case SPEX_CSC:
        {

            switch (A->kind)
            {

                //--------------------------------------------------------------
                // A is CSC, C is CSC
                //--------------------------------------------------------------

                case SPEX_CSC:
                {
                    // allocate C
                    SPEX_CHECK (SPEX_matrix_allocate (&C, SPEX_CSC, C_type,
                        m, n, nz, false, true, option)) ;
                    // copy the pattern of A into C
                    memcpy (C->p, A->p, (n+1) * sizeof (int64_t)) ;
                    memcpy (C->i, A->i, nz * sizeof (int64_t)) ;
                    // copy and typecast A->x into C->x
                    SPEX_CHECK (spex_cast_array (SPEX_X (C), C->type,
                        SPEX_X (A), A->type, nz, C->scale, A->scale, option)) ;
                }
                break ;

                //--------------------------------------------------------------
                // A is triplet, C is CSC
                //--------------------------------------------------------------

                case SPEX_TRIPLET:
                {

                    // Y = typecast the values of A into the type of C
                    // (not the pattern; Y is SPEX_DENSE)
                    SPEX_CHECK (spex_cast_matrix (&Y, C_type, A, option)) ;

                    // allocate workspace
                    W = (int64_t *) SPEX_calloc (n, sizeof (int64_t)) ;
                    if (W == NULL)
                    {
                        SPEX_FREE_ALL ;
                        return (SPEX_OUT_OF_MEMORY) ;
                    }

                    // allocate C
                    SPEX_CHECK (SPEX_matrix_allocate (&C, SPEX_CSC,
                        C_type, m, n, nz, false, true, option)) ;

                    // Scaling factor of C is currently in Y, set it
                    // here
                    SPEX_mpq_set(C->scale, Y->scale);
                    
                    // count the # of entries in each column
                    for (int64_t k = 0 ; k < nz ; k++)
                    {
                        W [A->j [k]]++ ;
                    }

                    // C->p = cumulative sum of W
                    spex_cumsum (C->p, W, n) ;

                    // build the matrix
                    switch (C->type)
                    {
                        case SPEX_MPZ:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t p = W [A->j [k]]++ ;
                                C->i [p] = A->i [k] ;
                                SPEX_CHECK (SPEX_mpz_set (
                                    SPEX_1D (C, p, mpz),
                                    SPEX_1D (Y, k, mpz))) ;
                            }
                            break ;

                        case SPEX_MPQ:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t p = W [A->j [k]]++ ;
                                C->i [p] = A->i [k] ;
                                SPEX_CHECK (SPEX_mpq_set (
                                    SPEX_1D (C, p, mpq),
                                    SPEX_1D (Y, k, mpq))) ;
                            }
                            break ;

                        case SPEX_MPFR:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t p = W [A->j [k]]++ ;
                                C->i [p] = A->i [k] ;
                                SPEX_CHECK (SPEX_mpfr_set (
                                    SPEX_1D (C, p, mpfr),
                                    SPEX_1D (Y, k, mpfr),
                                    round)) ;
                            }
                            break ;

                        case SPEX_INT64:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t p = W [A->j [k]]++ ;
                                C->i [p] = A->i [k] ;
                                SPEX_1D (C, p, int64) =
                                    SPEX_1D (Y, k, int64) ;
                            }
                            break ;

                        case SPEX_FP64:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t p = W [A->j [k]]++ ;
                                C->i [p] = A->i [k] ;
                                SPEX_1D (C, p, fp64) =
                                    SPEX_1D (Y, k, fp64) ;
                            }
                            break ;

                    }

                }
                break ;

                //--------------------------------------------------------------
                // A is dense, C is CSC
                //--------------------------------------------------------------

                case SPEX_DENSE:
                {
                    // Y = typecast the values of A into the type of C
                    SPEX_CHECK (spex_cast_matrix (&Y, C_type, A, option)) ;
                    int s ;

                    // count the actual nonzeros in Y
                    int64_t actual = 0 ;
                    switch (Y->type)
                    {

                        case SPEX_MPZ:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                SPEX_CHECK (SPEX_mpz_sgn (&s,
                                    SPEX_1D (Y, k, mpz))) ;
                                if (s != 0) actual++ ;
                            }
                            break ;

                        case SPEX_MPQ:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                SPEX_CHECK (SPEX_mpq_sgn (&s,
                                    SPEX_1D (Y, k, mpq))) ;
                                if (s != 0) actual++ ;
                            }
                            break ;

                        case SPEX_MPFR:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                SPEX_CHECK (SPEX_mpfr_sgn (&s,
                                    SPEX_1D (Y, k, mpfr))) ;
                                if (s != 0) actual++ ;
                            }
                            break ;

                        case SPEX_INT64:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                if (SPEX_1D (Y, k, int64) != 0) actual++ ;
                            }
                            break ;

                        case SPEX_FP64:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                if (SPEX_1D (Y, k, fp64) != 0) actual++ ;
                            }
                            break ;

                    }
                    // allocate C
                    SPEX_CHECK (SPEX_matrix_allocate (&C, SPEX_CSC, C_type,
                        m, n, actual, false, true, option)) ;
                        
                    // C's scaling factor is currently in Y. Set it here
                    SPEX_mpq_set(C->scale, Y->scale);

                    // Construct C
                    nz = 0 ;
                    switch (C->type)
                    {

                        case SPEX_MPZ:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                C->p [j] = nz ;
                                for (int64_t i = 0 ; i < m ; i++)
                                {
                                    SPEX_CHECK( SPEX_mpz_sgn( &s,
                                        Y->x.mpz[ i + j*A->m]));
                                    if (s != 0)
                                    {
                                        C->i [nz] = i ;
                                        SPEX_CHECK( SPEX_mpz_set (
                                            SPEX_1D (C, nz, mpz),
                                            Y->x.mpz[ i + j*A->m] ));
                                        nz++ ;
                                    }
                                }
                            }
                            break ;

                        case SPEX_MPQ:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                C->p [j] = nz ;
                                for (int64_t i = 0 ; i < m ; i++)
                                {
                                    SPEX_CHECK (SPEX_mpq_sgn (&s,
                                        Y->x.mpq[ i + j*A->m])) ;
                                    if (s != 0)
                                    {
                                        C->i [nz] = i ;
                                        SPEX_CHECK(SPEX_mpq_set (
                                            SPEX_1D(C, nz, mpq),
                                            Y->x.mpq[ i + j*A->m]));
                                        nz++ ;
                                    }
                                }
                            }
                            break ;

                        case SPEX_MPFR:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                C->p [j] = nz ;
                                for (int64_t i = 0 ; i < m ; i++)
                                {
                                    SPEX_CHECK (SPEX_mpfr_sgn (&s,
                                        Y->x.mpfr[i + j*A->m])) ;
                                    if (s != 0)
                                    {
                                        C->i [nz] = i ;
                                        SPEX_CHECK (SPEX_mpfr_set (
                                            SPEX_1D (C, nz, mpfr),
                                            Y->x.mpfr[i + j*A->m],
                                            round)) ;
                                        
                                        nz++ ;
                                    }
                                }
                            }
                            break ;

                        case SPEX_INT64:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                C->p [j] = nz ;
                                for (int64_t i = 0 ; i < m ; i++)
                                {
                                    if ( Y->x.int64[i +j*A->m] != 0)
                                    {
                                        C->i [nz] = i ;
                                        SPEX_1D (C, nz, int64) =
                                            Y->x.int64[i +j*A->m] ;
                                        nz++ ;
                                    }
                                }
                            }
                            break ;

                        case SPEX_FP64:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                C->p [j] = nz ;
                                for (int64_t i = 0 ; i < m ; i++)
                                {
                                    if ( Y->x.fp64[i +j*A->m] != 0)
                                    {
                                        C->i [nz] = i ;
                                        SPEX_1D (C, nz, fp64) =
                                            Y->x.fp64[i +j*A->m];
                                        nz++ ;
                                    }
                                }
                            }
                            break ;
                    }
                    C->p [n] = nz ;
                }
                break ;

                //--------------------------------------------------------------
                // A is dynamic_CSC, C is CSC
                //--------------------------------------------------------------

                case SPEX_DYNAMIC_CSC:
                {
                    // convert A to a SPEX_CSC x SPEX_MPZ matrix T
                    SPEX_CHECK(spex_dynamic_to_CSC_mpz(&T, A, nz, option));

                    if (C_type == SPEX_MPZ)
                    {
                        C = T;
                        T = NULL;
                    }
                    else
                    {
                        SPEX_CHECK (SPEX_matrix_copy (&C, SPEX_CSC, C_type,
                            T, option)) ;
                    }
                    SPEX_matrix_free (&T, option) ;
                }
                break;
            }

        }
        break ;

        //----------------------------------------------------------------------
        // C is triplet
        //----------------------------------------------------------------------

        case SPEX_TRIPLET:
        {

            switch (A->kind)
            {

                //--------------------------------------------------------------
                // A is CSC, C is triplet
                //--------------------------------------------------------------

                case SPEX_CSC:
                {
                    // allocate C
                    SPEX_CHECK (SPEX_matrix_allocate (&C, SPEX_TRIPLET, C_type,
                        m, n, nz, false, true, option)) ;
                    // copy and typecast A->x into C->x
                    SPEX_CHECK (spex_cast_array (SPEX_X (C), C->type,
                        SPEX_X (A), A->type, nz, C->scale, A->scale, option)) ;
                    // copy the row indices A->i into C->i
                    memcpy (C->i, A->i, nz * sizeof (int64_t)) ;
                    // construct C->j
                    for (int64_t j = 0 ; j < n ; j++)
                    {
                        for (int64_t p = A->p [j] ; p < A->p [j+1] ; p++)
                        {
                            C->j [p] = j ;
                        }
                    }
                    // set C->nz
                    C->nz = nz;
                }
                break ;

                //--------------------------------------------------------------
                // A is triplet, C is triplet
                //--------------------------------------------------------------

                case SPEX_TRIPLET:
                {
                    // allocate C
                    SPEX_CHECK (SPEX_matrix_allocate (&C, SPEX_TRIPLET, C_type,
                        m, n, nz, false, true, option)) ;
                    // copy the pattern of A into C
                    memcpy (C->j, A->j, nz * sizeof (int64_t)) ;
                    memcpy (C->i, A->i, nz * sizeof (int64_t)) ;
                    // copy and typecast A->x into C->x
                    SPEX_CHECK (spex_cast_array (SPEX_X (C), C->type,
                        SPEX_X (A), A->type, nz, C->scale, A->scale, option)) ;
                    // set C->nz
                    C->nz = nz;
                }
                break ;

                //--------------------------------------------------------------
                // A is dense, C is triplet
                //--------------------------------------------------------------

                case SPEX_DENSE:
                {
                    // convert A to a temporary CSC matrix
                    SPEX_CHECK (SPEX_matrix_copy (&T, SPEX_CSC, C_type,
                        A, option)) ;
                    // convert T from CSC to triplet
                    SPEX_CHECK (SPEX_matrix_copy (&C, SPEX_TRIPLET, C_type,
                        T, option)) ;
                    SPEX_matrix_free (&T, option) ;
                    // set C->nz
                    C->nz = nz;
                }
                break ;

                //--------------------------------------------------------------
                // A is dynamic_CSC, C is triplet
                //--------------------------------------------------------------

                case SPEX_DYNAMIC_CSC:
                {
                    // convert A to a SPEX_CSC x SPEX_MPZ matrix T
                    SPEX_CHECK(spex_dynamic_to_CSC_mpz(&T, A, nz, option));

                    SPEX_CHECK (SPEX_matrix_copy (&C, SPEX_TRIPLET, C_type,
                        T, option)) ;
                    SPEX_matrix_free (&T, option) ;
                }
                break;

            }

        }
        break ;

        //----------------------------------------------------------------------
        // C is dense
        //----------------------------------------------------------------------

        case SPEX_DENSE:
        {

            switch (A->kind)
            {

                //--------------------------------------------------------------
                // A is CSC, C is dense
                //--------------------------------------------------------------

                case SPEX_CSC:
                {
                    // allocate C
                    SPEX_CHECK (SPEX_matrix_allocate (&C, SPEX_DENSE, C_type,
                        m, n, nz, false, true, option)) ;

                    // Y = typecast the values of A into the type of C
                    SPEX_CHECK (spex_cast_matrix (&Y, C->type, A, option)) ;
                    
                    // Set C's scaling factor
                    SPEX_mpq_set(C->scale, Y->scale);

                    switch (C->type)
                    {

                        case SPEX_MPZ:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                for (int64_t p = A->p [j] ; p < A->p [j+1] ;p++)
                                {
                                    int64_t i = A->i [p] ;
                                    SPEX_CHECK (SPEX_mpz_set (
                                        SPEX_2D (C, i, j, mpz),
                                        SPEX_1D (Y, p, mpz))) ;
                                }
                            }
                            break ;

                        case SPEX_MPQ:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                for (int64_t p = A->p [j] ; p < A->p [j+1] ;p++)
                                {
                                    int64_t i = A->i [p] ;
                                    SPEX_CHECK (SPEX_mpq_set (
                                        SPEX_2D (C, i, j, mpq),
                                        SPEX_1D (Y, p, mpq))) ;
                                }
                            }
                            break ;

                        case SPEX_MPFR:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                for (int64_t p = A->p [j] ; p < A->p [j+1] ;p++)
                                {
                                    int64_t i = A->i [p] ;
                                    SPEX_CHECK (SPEX_mpfr_set (
                                        SPEX_2D (C, i, j, mpfr),
                                        SPEX_1D (Y, p, mpfr),
                                        round)) ;
                                }
                            }
                            break ;

                        case SPEX_INT64:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                for (int64_t p = A->p [j] ; p < A->p [j+1] ;p++)
                                {
                                    int64_t i = A->i [p] ;
                                    SPEX_2D (C, i, j, int64) =
                                        SPEX_1D (Y, p, int64) ;
                                }
                            }
                            break ;

                        case SPEX_FP64:
                            for (int64_t j = 0 ; j < n ; j++)
                            {
                                for (int64_t p = A->p [j] ; p < A->p [j+1] ;p++)
                                {
                                    int64_t i = A->i [p] ;
                                    SPEX_2D (C, i, j, fp64) =
                                        SPEX_1D (Y, p, fp64) ;
                                }
                            }
                            break ;

                    }

                }
                break ;

                //--------------------------------------------------------------
                // A is triplet, C is dense
                //--------------------------------------------------------------

                case SPEX_TRIPLET:
                {
                    // allocate C
                    SPEX_CHECK (SPEX_matrix_allocate (&C, SPEX_DENSE, C_type,
                        m, n, nz, false, true, option)) ;

                    // Y = typecast the values of A into the type of C
                    SPEX_CHECK (spex_cast_matrix (&Y, C->type, A, option)) ;

                    // Set C's scaling factor
                    SPEX_mpq_set(C->scale, Y->scale);
                    switch (C->type)
                    {

                        case SPEX_MPZ:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t i = A->i [k] ;
                                int64_t j = A->j [k] ;
                                SPEX_CHECK (SPEX_mpz_set (
                                    SPEX_2D (C, i, j, mpz),
                                    SPEX_1D (Y, k, mpz))) ;
                            }
                            break ;

                        case SPEX_MPQ:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t i = A->i [k] ;
                                int64_t j = A->j [k] ;
                                SPEX_CHECK (SPEX_mpq_set (
                                    SPEX_2D (C, i, j, mpq),
                                    SPEX_1D (Y, k, mpq))) ;
                            }
                            break ;

                        case SPEX_MPFR:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t i = A->i [k] ;
                                int64_t j = A->j [k] ;
                                SPEX_CHECK (SPEX_mpfr_set (
                                    SPEX_2D (C, i, j, mpfr),
                                    SPEX_1D (Y, k, mpfr),
                                    round)) ;
                            }
                            break ;

                        case SPEX_INT64:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t i = A->i [k] ;
                                int64_t j = A->j [k] ;
                                SPEX_2D (C, i, j, int64) =
                                    SPEX_1D (Y, k, int64) ;
                            }
                            break ;

                        case SPEX_FP64:
                            for (int64_t k = 0 ; k < nz ; k++)
                            {
                                int64_t i = A->i [k] ;
                                int64_t j = A->j [k] ;
                                SPEX_2D (C, i, j, fp64) =
                                    SPEX_1D (Y, k, fp64) ;
                            }
                            break ;

                    }
                }
                break ;

                //--------------------------------------------------------------
                // A is dense, C is dense
                //--------------------------------------------------------------

                case SPEX_DENSE:
                {
                    // allocate C
                    SPEX_CHECK (SPEX_matrix_allocate (&C, SPEX_DENSE, C_type,
                        m, n, nz, false, true, option)) ;

                    // copy and typecast A->x into C->x
                    SPEX_CHECK (spex_cast_array (SPEX_X (C), C->type,
                        SPEX_X (A), A->type, nz, C->scale, A->scale, option)) ;
                }
                break ;

                //--------------------------------------------------------------
                // A is dynamic_CSC, C is dense
                //--------------------------------------------------------------

                case SPEX_DYNAMIC_CSC:
                {
                    // convert A to a SPEX_CSC x SPEX_MPZ matrix T
                    SPEX_CHECK(spex_dynamic_to_CSC_mpz(&T, A, nz, option));

                    SPEX_CHECK (SPEX_matrix_copy (&C, SPEX_DENSE, C_type,
                        T, option)) ;
                    SPEX_matrix_free (&T, option) ;
                }
                break;

            }

        }
        break ;

        //----------------------------------------------------------------------
        // C is dynamic_CSC
        //----------------------------------------------------------------------

        case SPEX_DYNAMIC_CSC:
        {
            if (A->kind != SPEX_DYNAMIC_CSC)
            {
                if (A->kind == SPEX_CSC && A->type == SPEX_MPZ)
                {
                    // A is already a CSC x MPZ matrix
                    SPEX_CHECK(spex_CSC_mpz_to_dynamic(&C, A, option));
                }
                else
                {
                    // convert A to a SPEX_CSC x SPEX_MPZ matrix T
                    SPEX_CHECK (SPEX_matrix_copy (&T, SPEX_CSC, SPEX_MPZ,
                        A, option)) ;
                    SPEX_CHECK(spex_CSC_mpz_to_dynamic(&C, T, option));
                    SPEX_matrix_free (&T, option) ;
                }
            }
            else // make a exact same copy
            {
                // allocate space for C
                SPEX_CHECK (SPEX_matrix_allocate (&C, SPEX_DYNAMIC_CSC,
                    C_type, m, n, 0, false, true, option)) ;
                // copy A->x into C->x
                for (int64_t j = 0; j < n; j++)
                {
                    int64_t Aj_nz = A->v[j]->nz;
                    // reallocate space for each column of C
                    SPEX_CHECK(SPEX_vector_realloc(C->v[j], Aj_nz, option));
                    memcpy (C->v[j]->i, A->v[j]->i, Aj_nz*sizeof (int64_t));

                    for (int64_t p = 0; p < Aj_nz; p++)
                    {
                        SPEX_CHECK(SPEX_mpz_set(C->v[j]->x[p],
                            A->v[j]->x[p]));
                    }
                    C->v[j]->nz = Aj_nz;
                }
            }

        }
        break ;

    }

    //--------------------------------------------------------------------------
    // free workspace and return result
    //--------------------------------------------------------------------------

    SPEX_FREE_WORK ;
    (*C_handle) = C ;

    return (SPEX_OK) ;
}
