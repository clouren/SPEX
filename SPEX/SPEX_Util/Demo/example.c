//------------------------------------------------------------------------------
// SPEX_Util/Demo/example.c: example for creating and use safely wrapped gmp
// functions
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


#include "my_safe_gmp.h"


/* This example shows how to extend the wrapper mechanism provided by SPEX to
 * additional gmp/mpfr functions and thus use them safely functions.
 */

// usage:
// example > out

#define FREE_WORKSPACE                       \
{                                            \
    my_safe_mpz_clear(z);                    \
    my_safe_mpq_clear(q);                    \
    SPEX_finalize() ;                        \
}

#define OK(method)                      \
{                                       \
    ok = method ;                       \
    if (ok != SPEX_OK)                  \
    {                                   \
        printf ("Error: %d line %d file %s\n", ok, __LINE__, __FILE__) ; \
        FREE_WORKSPACE ;                \
        return 0 ;                      \
    }                                   \
}

   
    
int main (void)
{
    
    //--------------------------------------------------------------------------
    // Prior to using SPEX, its environment must be initialized. This is done
    // by calling the SPEX_initialize() function.
    //--------------------------------------------------------------------------

    SPEX_initialize();

    //--------------------------------------------------------------------------
    // Declare and initialize essential variables
    //--------------------------------------------------------------------------

    printf ("\nThis example shows how to extend the wrapper mechanism ");
    printf ("provided by SPEX to additional gmp/mpfr functions and thus ");
    printf ("use them safely functions.\n\n") ;
    SPEX_info ok;
    mpz_t z;
    mpq_t q;
    if (SPEX_mpz_init(z) != SPEX_OK)
    {
        SPEX_finalize() ;
        return 0;
    }
    if (SPEX_mpq_init(q) != SPEX_OK)
    {
        my_safe_mpz_clear(z);
        SPEX_finalize() ;
        return 0;
    }

    //--------------------------------------------------------------------------
    // safely use some gmp/mpfr functions
    //--------------------------------------------------------------------------
    // q = 1/2
    OK(SPEX_mpq_set_ui(q, 1, 2));
    // z = den(q) = 2
    OK(my_safe_mpq_get_den(z, q));

    // print z and q
    OK(my_safe_gmp_printf("q = %Qd\nz = %Zd\n", q, z));

    //--------------------------------------------------------------------------
    // Free memory
    //--------------------------------------------------------------------------

    FREE_WORKSPACE;
    printf ("\n%s: all tests passed\n\n", __FILE__) ;
    return 0;
}

