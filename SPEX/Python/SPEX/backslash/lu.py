# SPEX: (c) 2022, Chris Lourenco, United States Naval Academy,
# Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
# Texas A&M University. All Rights Reserved.
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

from SPEX.utils import Options
from SPEX._connect import spex_connect

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import scipy
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix, isspmatrix, isspmatrix_csc, linalg

def lu( A, b, options=Options('double', 'colamd')):
    ## A is a scipy.sparse(data must be float64) #technically it only needs to be numerical
    ## b is a numpy.array (data must be float64)
    ## options is a dictionary that specifies what tipe the solution should be, this by default is double

    ##--------------------------------------------------------------------------
    ## Verify inputs
    ##--------------------------------------------------------------------------
    if not isspmatrix(A):
        print("Input matrix must be scipy.sparse")
        raise TypeError
    ## If the sparse input matrix is not in csc form, convert it into csc form
    if not isspmatrix_csc(A):
        A.tocsc()
    # Check input shape
    if A.shape[1]!=b.shape[0]:
        print("Matrix input shapes must match")
        raise TypeError

    if options.ordering==None:
        options.default_lu()

    ##--------------------------------------------------------------------------
    ## Call SPEX
    ##--------------------------------------------------------------------------
    x=spex_connect(A,b,options.order(),options.charOut(),2) #3 calls lu

    return x

