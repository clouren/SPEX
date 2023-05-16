//------------------------------------------------------------------------------
// SPEX_QR/Source/spex_qr_ipgs.c: Integer Preserving Gram-Schmidt
//------------------------------------------------------------------------------

// SPEX_QR: (c) 2021-2023, Chris Lourenco, Lorena Mejia Domenzain,
// Timothy A. Davis, and Erick Moreno-Centeno. All Rights Reserved.
// SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

//------------------------------------------------------------------------------


/* This code performs one iteration of REF QR via Integer-preserving Gram-Schmidt
 */

# include "spex_qr_internal.h"


SPEX_info SPEX_qr_factorize
(
    SPEX_matrix *R_handle,    // Null on input, contains R on output
    SPEX_matrix *Q_handle,    // Null on input, contains Q on output
    SPEX_matrix rhos,         // sequence of pivots
    const SPEX_matrix A,      // Matrix to be factored
    SPEX_options option
)
{

    //Copy A into Q
    for(i=0;i<A->nnz;i++)
    {
        
    }

    // Perform IPGS to get Q and R
    for (k=0;k<n;k++)
    {
        SPEX_CHECK(spex_qr_ipge());
    }

}