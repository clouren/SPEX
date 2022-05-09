//------------------------------------------------------------------------------
// SPEX_Util/Demo/my_safe_gmp.c: wrapped gmp functions for the demo programs
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


#include "my_safe_gmp.h"

// ignore warnings about unused parameters in this file
#pragma GCC diagnostic ignored "-Wunused-parameter"

//------------------------------------------------------------------------------
// my_safe_gmp_printf
//------------------------------------------------------------------------------

/* Safely print to the standard output stdout. Return positive value (the number
 * of characters written) upon success, otherwise return negative value (error
 * code) */

SPEX_info my_safe_gmp_printf
(
    const char *format,
    ...
)
{
    // use this macro since no memory space will be allocated
    // Start the GMP wrapper
    SPEX_GMP_WRAPPER_START ;
    
    // call gmp_vprintf 
    va_list args; 
    va_start (args, format) ;
    int n = gmp_vprintf (format, args) ; 
    va_end (args) ;

    // Finish the wrapper
    SPEX_GMP_WRAPPER_FINISH ;
    // gmp_vprintf returns -1 if an error occurred.
    return ((n < 0) ? SPEX_INCORRECT_INPUT : SPEX_OK) ;
}

//------------------------------------------------------------------------------
// my_safe_mpq_get_den
//------------------------------------------------------------------------------

/* Purpose: Safely set an mpz number = denominator of an mpq number */

SPEX_info my_safe_mpq_get_den
(
    mpz_t x,
    const mpq_t y
)
{
    // use this macro to save the address of x. In case of failure, and newly
    // allocated memory space for x will be free'd, so caller of this function
    // need not to aware of that.
    SPEX_GMPZ_WRAPPER_START (x) ;
    mpq_get_den (x, y) ;
    SPEX_GMP_WRAPPER_FINISH ;
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// my_safe_mpz_clear
//------------------------------------------------------------------------------

/* Purpose: Safely free memory space allocated for a mpz object */

SPEX_info my_safe_mpz_clear
(
    mpz_t x
)
{
    // use this macro since no memory space will be allocated
    SPEX_GMP_WRAPPER_START;
    mpz_clear (x) ;
    SPEX_GMP_WRAPPER_FINISH ;
    return (SPEX_OK) ;
}

//------------------------------------------------------------------------------
// my_safe_mpq_clear
//------------------------------------------------------------------------------

/* Purpose: Safely free memory space allocated for a mpq object */

SPEX_info my_safe_mpq_clear
(
    mpq_t x
)
{
    // use this macro since no memory space will be allocated
    SPEX_GMP_WRAPPER_START ;
    mpq_clear (x) ;
    SPEX_GMP_WRAPPER_FINISH ;
    return (SPEX_OK) ;
}


