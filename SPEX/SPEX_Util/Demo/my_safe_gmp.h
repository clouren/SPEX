//------------------------------------------------------------------------------
// SPEX_Util/Demo/my_safe_gmp.h: #include file the demo programs
//------------------------------------------------------------------------------

// SPEX_Util: (c) 2019-2021, Chris Lourenco (US Naval Academy), Jinhao Chen,
// Erick Moreno-Centeno, Timothy A. Davis, Texas A&M.  All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------
#ifndef MY_SAFE_GMP_H
#define MY_SAFE_GMP_H

#include "SPEX.h"
#include "spex_gmp.h"

/* Safely print to the standard output stdout. Return positive value (the number
 * of characters written) upon success, otherwise return negative value (error
 * code) */
SPEX_info my_safe_gmp_printf (const char *format, ... ) ;

/* Purpose: Safely set an mpz number = denominator of an mpq number */
SPEX_info my_safe_mpq_get_den (mpz_t x, const mpq_t y) ;

/* Purpose: Safely free memory space allocated for a mpz object */
SPEX_info my_safe_mpz_clear(mpz_t x);

/* Purpose: Safely free memory space allocated for a mpq object */
SPEX_info my_safe_mpq_clear(mpq_t x);

#endif
