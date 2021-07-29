import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix, isspmatrix, isspmatrix_csc
#TODO add options
#TODO add the functionality of getting strings as output (related to options, so user can specify)
#TODO each function in a different file? (check how that works in python libraries)


def SPEX_Chol_backslashVIEJO( A,b ): 
    ## A is a scipy.sparse.csc_matrix (data must be float64) #technically it only needs to be numerical
    ## b is a numpy.array (data must be float64)
    
    ##--------------------------------------------------------------------------
    ## Load the library with the "C bridge code"
    ##--------------------------------------------------------------------------
    lib = ctypes.CDLL('./libraryE.so')
    c_backslash = lib.python_backslash2
    
    ##--------------------------------------------------------------------------
    ## Specify the parameter types and return type of the C function
    ##--------------------------------------------------------------------------
    c_backslash.argtypes = [ctypes.POINTER(ctypes.c_double),
                            ndpointer(dtype=np.int64, ndim=1, flags=None),
                            ndpointer(dtype=np.int64, ndim=1, flags=None),
                            ndpointer(dtype=np.float64, ndim=1, flags=None),
                            ndpointer(dtype=np.float64, ndim=1, flags=None), 
                            ctypes.c_int, 
                            ctypes.c_int]
    c_backslash.restype = None #C method is void
    
    n=A.shape[0] #number of columns/rows of A
    
    x = (ctypes.c_double*n)() #empty solution vector

    ##--------------------------------------------------------------------------
    ## Solve Ax=b using REF Sparse Cholesky Factorization
    ##--------------------------------------------------------------------------
    c_backslash(x,  
                A.indptr.astype(np.int64), #without the cast it would be int32 and it would not be compatible with the C method
                A.indices.astype(np.int64),
                A.data.astype(np.float64), 
                b,
                n,
                A.nnz)
   
    return x

def SPEX_Chol_string( A, b, order ): 
    ## A is a scipy.sparse.csc_matrix (data must be float64) #technically it only needs to be numerical
    ## b is a numpy.array (data must be float64)
    
    ##--------------------------------------------------------------------------
    ## Load the library with the "C bridge code"
    ##--------------------------------------------------------------------------
    lib = ctypes.CDLL('./SPEX_Chol_connect.so')
    c_backslash = lib.python_backslash_char
    
    ##--------------------------------------------------------------------------
    ## Specify the parameter types and return type of the C function
    ##--------------------------------------------------------------------------
    c_backslash.argtypes = [ctypes.POINTER(ctypes.POINTER(ctypes.c_char)),
                            ndpointer(dtype=np.int64, ndim=1, flags=None),
                            ndpointer(dtype=np.int64, ndim=1, flags=None),
                            ndpointer(dtype=np.float64, ndim=1, flags=None),
                            ndpointer(dtype=np.float64, ndim=1, flags=None), 
                            ctypes.c_int, 
                            ctypes.c_int]
    c_backslash.restype = None #C method is void
    
    n=A.shape[0] #number of columns/rows of A
    
    x = ctypes.POINTER(ctypes.c_char)() #empty solution vector
    print(type(x))

    ##--------------------------------------------------------------------------
    ## Solve Ax=b using REF Sparse Cholesky Factorization
    ##--------------------------------------------------------------------------
    c_backslash(x,  
                A.indptr.astype(np.int64), #without the cast it would be int32 and it would not be compatible with the C method
                A.indices.astype(np.int64),
                A.data.astype(np.float64), 
                b,
                n,
                A.nnz)
   
    print(type(x))
    return x

def SPEX_Chol_double( A, b, order ): 
    ## A is a scipy.sparse.csc_matrix (data must be float64) #technically it only needs to be numerical
    ## b is a numpy.array (data must be float64)
    
    ##--------------------------------------------------------------------------
    ## Load the library with the "C bridge code"
    ##--------------------------------------------------------------------------
    lib = ctypes.CDLL('./SPEX_Chol_connect.so')
    c_backslash = lib.python_backslash_double
    
    ##--------------------------------------------------------------------------
    ## Specify the parameter types and return type of the C function
    ##--------------------------------------------------------------------------
    c_backslash.argtypes = [ctypes.POINTER(ctypes.c_double),
                            ndpointer(dtype=np.int64, ndim=1, flags=None),
                            ndpointer(dtype=np.int64, ndim=1, flags=None),
                            ndpointer(dtype=np.float64, ndim=1, flags=None),
                            ndpointer(dtype=np.float64, ndim=1, flags=None), 
                            ctypes.c_int, 
                            ctypes.c_int]
    c_backslash.restype = None #C method is void
    
    n=A.shape[0] #number of columns/rows of A
    
    x = (ctypes.c_double*n)() #empty solution vector

    ##--------------------------------------------------------------------------
    ## Solve Ax=b using REF Sparse Cholesky Factorization
    ##--------------------------------------------------------------------------
    c_backslash(x,  
                A.indptr.astype(np.int64), #without the cast it would be int32 and it would not be compatible with the C method
                A.indices.astype(np.int64),
                A.data.astype(np.float64), 
                b,
                n,
                A.nnz)
   
    return np.array(list(x), dtype=np.float64)

def Cholesky( A, b, options={'SolutionType': 'double', 'Ordering': 'amd'}): 
    ## A is a scipy.sparse(data must be float64) #technically it only needs to be numerical
    ## b is a numpy.array (data must be float64)
    ## options is a dictionary that specifies what tipe the solution should be, this by default is double
    
    ##--------------------------------------------------------------------------
    ## Verify inputs
    ##--------------------------------------------------------------------------
    if not isspmatrix(A):
        return "Error: A must be a scipy sparse matrix"
    ## If the sparse input matrix is not in csc form, convert it into csc form
    if not isspmatrix_csc(A):
        A.tocsc()
    
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
        return "Error: invalid options"
        
    ##--------------------------------------------------------------------------
    ## Call the correct function depending on the desired output type
    ##--------------------------------------------------------------------------    
    if options['SolutionType']=="double":
        x=SPEX_Chol_double(A,b,order)
    elif options['SolutionType']=="string":
        x=SPEX_Chol_string(A,b,order)
    else:
        return "Error: invalid options"

   
    return x


def spex_matrix_from_file(fname):
    #fname is the name of the file that contains matrix A
  
    ##--------------------------------------------------------------------------
    ## Load file data and store it in mat
    ##--------------------------------------------------------------------------
    mat = np.loadtxt(fname, delimiter=' ')
    
    ##--------------------------------------------------------------------------
    ## Splice mat to get the componets of matrix A and populate a scipy sparse matrix
    ##--------------------------------------------------------------------------
    data = mat[1:,2]
    row = mat[1:,0].astype(int)
    col = mat[1:,1].astype(int)
    n = mat[0][0].astype(int)
    m = mat[0][1].astype(int)
    ## The matrix in the file is in triplet form
    triplet=coo_matrix((data, (row, col)),shape=(n, m))
    ## Change the triplet matrix into a csc matrix
    A = triplet.tocsc()
    
    return A


