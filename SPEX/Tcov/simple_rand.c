// ----------------------------------------------------------------------------
// SPEX/Tcov/simple_rand.c: a very simple random number generator
// ----------------------------------------------------------------------------

// SPEX: (c) 2019-2024, Christopher Lourenco, Jinhao Chen,
// Lorena Mejia Domenzain, Erick Moreno-Centeno, and Timothy A. Davis.
// All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//-----------------------------------------------------------------------------

//  The POSIX.1-2001 example of rand, duplicated here so that the same sequence
//  will be generated on different machines.  The purpose is not to generate
//  high-quality random numbers, but to ensure the same sequence is generated
//  on any computer or operating system.  This allows the test coverage to be
//  repeatable.

//  Since the simple_rand ( ) function is replicated from the POSIX.1-2001
//  standard, no copyright claim is intended for this specific file.  The
//  copyright statement above applies to all of SPEX, not this file.

#include "simple_rand.h"

/* simple_rand is not thread-safe */
uint64_t simple_rand_next = 1 ;

#define SIMPLE_RAND_MAX 32767

/* return a random number between 0 and SIMPLE_RAND_MAX */
uint64_t simple_rand (void)
{
   simple_rand_next = simple_rand_next * 1103515245 + 12345 ;
   return ((simple_rand_next/65536) % (SIMPLE_RAND_MAX + 1));
}

/* set the seed */
void simple_rand_seed (uint64_t seed)
{
   simple_rand_next = seed ;
}

/* get the seed */
uint64_t simple_rand_getseed (void)
{
   return (simple_rand_next);
}

/* return a random double between 0 and 1, inclusive */
double simple_rand_x (void )
{
    return (((double) simple_rand ( )) / ((double) SIMPLE_RAND_MAX));
}

/* return a random uint64_t */
uint64_t simple_rand_i (void )
{
    uint64_t i = 0 ;
    for (int k = 0 ; k < 5 ; k++)
    {
        i = SIMPLE_RAND_MAX * i + simple_rand ( );
    }
    return (i);
}

