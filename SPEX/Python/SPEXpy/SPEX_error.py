#-------------------------------------------------------------------------------
# SPEX/Python/SPEXpy/SPEX_error.py: class SPEX_error
#-------------------------------------------------------------------------------

# SPEX: (c) 2022, Chris Lourenco, Jinhao Chen,
# Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
# All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

#------------------------------------------------------------------------------

class SPEX_error(LookupError):
    '''raise this when there's a lookup error for spex'''


def determine_error(ok):
    errorMessages={
        1:"out of memory",
        2:"the input matrix A is singular",
        3:"one or more input arguments are incorrect",
        4:"the input matrix is unsymmetric",
        5:"the input matrix is not SPD",
        6:"the algorithm is not compatible with the factorization",
        7:"SPEX used without proper initialization",
    }
    return errorMessages.get(ok*(-1))
