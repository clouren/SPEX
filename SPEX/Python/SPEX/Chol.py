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

def backslash( A, b, options={'SolutionType': 'double', 'Ordering': 'amd'}): 
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
    ## Check symmetry    
    tol=1e-8    
    if scipy.sparse.linalg.norm(A-A.T, scipy.Inf) > tol:
        print("Input matrix is not symmetric")
        raise TypeError
    
    ##--------------------------------------------------------------------------
    ## Ordering
    ##--------------------------------------------------------------------------
    if options['Ordering']=="none":
        order=0
    elif options['Ordering']=="colamd":
        order=1
    elif options['Ordering']=="amd": ##amd is the default ordering for Cholesky
        order=2
    else:
        print("Invalid order options")
        raise ValueError
        
    ##--------------------------------------------------------------------------
    ## Call the correct function depending on the desired output type
    ##--------------------------------------------------------------------------    
    if options['SolutionType']=="double":
        charOut=0
    elif options['SolutionType']=="string":
        charOut=0
    else:
        print("Invalid output type options")
        raise ValueError

    x=spex_chol_backslash(A,b,order,charOut)

    return x

def spex_chol_backslash( A, b, order, charOut ): 
    ## A is a scipy.sparse.csc_matrix (data must be float64) #technically it only needs to be numerical
    ## b is a numpy.array (data must be float64)
    
    ##--------------------------------------------------------------------------
    ## Load the library with the "C bridge code"
    ##--------------------------------------------------------------------------
    lib = ctypes.CDLL('./SPEX_connect.so')
    c_backslash = lib.spex_chol_python
    
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
                            ctypes.c_bool]
    c_backslash.restype = ctypes.c_int #C method is void
    
    n=A.shape[0] #number of columns/rows of A
    
    x_v = (ctypes.c_void_p*n)()
        
    ##--------------------------------------------------------------------------
    ## Solve Ax=b using REF Sparse Cholesky Factorization
    ##--------------------------------------------------------------------------
    c_backslash(x_v,
                A.indptr.astype(np.int64), #without the cast it would be int32 and it would not be compatible with the C method
                A.indices.astype(np.int64),
                A.data.astype(np.float64), 
                b,
                n,
                n,
                A.nnz,
                order,
                charOut)
    
    ##--------------------------------------------------------------------------
    ## Cast solution into correct type (string or double)
    ##--------------------------------------------------------------------------
    if charOut:
        x = ctypes.cast(x_v, ctypes.POINTER(ctypes.c_char_p))
    else:
        #x = ctypes.cast(x_v, ctypes.POINTER(ctypes.c_double))
        x=[]
        for i in range(n):
            val=ctypes.cast(x_v[i], ctypes.POINTER(ctypes.c_double))
            x.append(val[0]) ##this can also be changed to be a numpy array instead of a list
    
    return x
