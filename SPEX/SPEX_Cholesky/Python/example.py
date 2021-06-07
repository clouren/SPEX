import ctypes
import numpy as np
from numpy.ctypeslib import ndpointer
from scipy.sparse import csc_matrix

#TODO add options
#TODO add the functionality of getting strings as output (related to options, so user can specify)

def SPEX_Chol_backslash( A,b ): 
    ## A is a scipy.sparse.csc_matrix (data must be float64) #technically it only needs to be numerical
    ## b is a numpy.array (data must be float64)
    
    # "connect" to the C function
    lib = ctypes.CDLL('./libraryE.so')
    c_backslash = lib.python_backslash2
    
    ## Specify the parameter types and return type of the C function
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

    ## Solve Ax=b using REF Sparse Cholesky Factorization
    c_backslash(x,  
                A.indptr.astype(np.int64), #without the cast it would be int32 and it would not be compatible with the C method
                A.indices.astype(np.int64),
                A.data.astype(np.float64), 
                b,
                n,
                A.nnz)
   
    return x


def test():
    row = np.array([0, 0, 0, 1, 1, 1, 2, 2, 2])
    col = np.array([0, 1, 2, 0, 1, 2, 0, 1, 2])
    data = np.array([4, 12, -16, 12, 37, -43, -16, -43, 98],dtype=np.float64)
    A=csc_matrix((data, (row, col)), shape=(3, 3))
    b=np.ones(3,dtype=np.float64)

    x=SPEX_Chol_backslash(A,b)
    print(list(x))

test()


def spex_matrix_from_file(fname):
    # Setting name of the file that the data is to be extarcted from in python
  
    # Loading file data into numpy array and storing it in variable called data_collected
    mat = np.loadtxt(fname, delimiter=' ')
    data = mat[1:,2]
    row = mat[1:,0].astype(int)
    col = mat[1:,1].astype(int)
    n = mat[0][0].astype(int)+1
    m = mat[0][1].astype(int)+1
    triplet=coo_matrix((data, (row, col)),shape=(n, m))
    
    csc = triplet.tocsc()
    
    return csc

def test2():
    fname='../ExampleMats/2.mat.txt
    A=spex_matrix_from_file(fname)
    b=np.ones(n,dtype=np.float64)
    x=SPEX_Chol_backslash.SPEX_Chol_backslash(A,b)
    print('yay')
    
test2()




