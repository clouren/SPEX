SPEX is a software package for SParse EXact algebra

Files and folders in this distribution:

    README.md           this file
    
    AMD                 Approximate minimum degree ordering for
                        preordering matrices prior to factorization.
                        Please see www.suitesparse.com for more information
    
    COLAMD              Column approximate minimum degree ordering
                        for preordering matrices prior to factorization.
                        Please see www.suitesparse.com for more information
                    
    include             Header files for SPEX that are generated by
                        the general makefile
    
    share               Share files for SPEX that are generated by the 
                        general makefile
    
    SPEX                Folder containing the SPEX factorization routines
                    
    SuiteSparse_config  Configure file for all SuiteSparse functions.
                        Please see www.suitesparse.com for information
    
    SPEX_UTIL   Utility functions for all SPEX components
    
    Makefile    compiles SPEX and its dependencies

Dependencies of SPEX:

    AMD                 approximate minimum degree ordering
    
    COLAMD              column approximate minimum degree ordering
    
    SuiteSparse_config  configuration for all of SuiteSparse
    
    GNU GMP             GNU Multiple Precision Arithmetic Library 
                        for big integer operations
    
    GNU MPFR            GNU Multiple Precision Floating-Point Reliable
                        Library for arbitrary precision floating point
                        operations

Default instalation locations:

    include
    lib
    share
    
To compile SPEX and its dependencies, just type cd to the SPEX folder and 
type "make"
This will also run a few short demos
To install the package system-wide, copy the `lib/*` to /usr/local/lib,
and copy `include/*` to /usr/local/include.

Primary Author: Chris Lourenco

Coauthors (alphabetical order):

    Jinhao Chen
    Tim Davis    
    Lorena Mejia Domenzain
    Erick Moreno-Centeno

