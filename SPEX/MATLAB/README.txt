TODO: Create a general SPEX matlab function.

x = SPEX_backslash(A,b)

if (unsymmetric): x = SPEX_Left_LU_backslash (A, b)

if (SPD): x = SPEX_Chol(A,b)

... more to come

