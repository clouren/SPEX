import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
from scipy.sparse import csc_matrix
from scipy.sparse import coo_matrix

#TODO add options
#TODO add the functionality of getting strings as output (related to options, so user can specify)
##TODO figure out paths
##TODO change file names :D
# TODO make output a np.array ...maybe
#TODO each function in a different file? (check how that works in python libraries)

##TODO add ifs so that A can be just a numpy matrix and then we make it csc

def SPEX_Chol_backslash( A,b ): 
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


def test():
    ##baby test with small dense matrix
    ##just to check the "connection" with c works
    row = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
    col = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
    data = np.array([4, 12, -16, 12, 37, -43, -16, -43, 98],dtype=np.float64)
    A=csc_matrix((data, (row, col)), shape=(3, 3))
    b=np.ones(3,dtype=np.float64)

    x=SPEX_Chol_backslash(A,b)
    print(list(x))
    xArray=np.array(list(x))
    print(xArray)
    


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


