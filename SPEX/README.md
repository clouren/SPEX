SPEX is a software package for SParse EXact algebra

Files and folders in this distribution:

    README.md   this file
    SPEX_Cholesky   Sparse integer-preserving SPEX_Cholesky
                    factorization for exactly solving SPD
                    linear systems
    SPEX_Left_LU    Sparse left-looking integer-preserving
                    LU factorization for exactly solve 
                    sparse linear systems
    SPEX_UTIL   Utility functions for all SPEX components
    Makefile    compiles SPEX and its dependencies

Dependencies (all part of SuiteSparse):

    AMD                 approximate minimum degree ordering
    COLAMD              column approximate minimum degree ordering
    SuiteSparse_config  configuration for all of SuiteSparse

Default instalation locations:

    include
    lib
    share
    
To compile SPEX and its dependencies, just type "make" in this folder.
This will also run a few short demos
To install the package system-wide, copy the `lib/*` to /usr/local/lib,
and copy `include/*` to /usr/local/include.

Primary Author: Chris Lourenco

Coauthors (alphabetical order):

    Jinhao Chen
    Tim Davis
    Erick Moreno-Centeno

