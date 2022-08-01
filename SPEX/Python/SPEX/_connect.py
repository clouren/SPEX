# SPEX: (c) 2022, Chris Lourenco, United States Naval Academy, 
# Lorena Mejia Domenzain, Jinhao Chen, Erick Moreno-Centeno, Timothy A. Davis,
# Texas A&M University. All Rights Reserved. 
# SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later

import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
import scipy
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix, isspmatrix, isspmatrix_csc, linalg

def spex_connect( A, b, order, charOut, algorithm ): 
    ## A is a scipy.sparse.csc_matrix (data must be float64) #technically it only needs to be numerical
    ## b is a numpy.array (data must be float64)
    
    ##--------------------------------------------------------------------------
    ## Load the library with the "C bridge code"
    ##--------------------------------------------------------------------------
    lib = ctypes.CDLL('./SPEX_connect.so')
    c_backslash = lib.spex_python
    
    ##--------------------------------------------------------------------------
    ## Specify the parameter types and return type of the C function
    ##--------------------------------------------------------------------------
    c_backslash.argtypes = [ctypes.POINTER(ctypes.c_void_p),
                            ndpointer(dtype=np.int64, ndim=1, flags=None),
                            ndpointer(dtype=np.int64, ndim=1, flags=None),
                            ndpointer(dtype=np.float64, ndim=1, flags=None),
                            ndpointer(dtype=np.float64, ndim=1, flags=None), 
                            ctypes.c_int, 
                            ctypes.c_int,
                            ctypes.c_int,
                            ctypes.c_int,
                            ctypes.c_int,
                            ctypes.c_bool]
    c_backslash.restype = ctypes.c_int
    
    n=A.shape[0] #number of columns/rows of A
    
    x_v = (ctypes.c_void_p*n)()
        
    ##--------------------------------------------------------------------------
    ## Solve Ax=b using REF Sparse Cholesky Factorization
    ##--------------------------------------------------------------------------
    ok=c_backslash(x_v,
                A.indptr.astype(np.int64), #without the cast it would be int32 and it would not be compatible with the C method
                A.indices.astype(np.int64),
                A.data.astype(np.float64), 
                b,
                n,
                n,
                A.nnz,
                order,
                algorithm,
                charOut)
    
    if ok!=0: 
        raise SPEXerror(determine_error(ok))
        
    ##--------------------------------------------------------------------------
    ## Cast solution into correct type (string or double)
    ##--------------------------------------------------------------------------
    if charOut:
        x = ctypes.cast(x_v, ctypes.POINTER(ctypes.c_char_p))
        x = castSol(x,n)
    else:
        #x = ctypes.cast(x_v, ctypes.POINTER(ctypes.c_double))
        x=[]
        for i in range(n):
            val=ctypes.cast(x_v[i], ctypes.POINTER(ctypes.c_double))
            x.append(val[0]) ##this can also be changed to be a numpy array instead of a list
    
    return np.array(x)

def castSol(val,n):
    x=[]
    for i in range(n):
        x.append(val[i])
    return x   

class SPEXerror(LookupError):
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