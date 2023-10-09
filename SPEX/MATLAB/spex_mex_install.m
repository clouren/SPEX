function spex_mex_install(run_demo)
% spex_mex_INSTALL: install and test the MATLAB interface to SPEX MATLAB functions.
%
% Usage: spex_mex_install
%
% Required Libraries: GMP, MPFR, AMD, COLAMD, SPEX.  If -lamd and -lcolamd are
% not available, install them with 'make install' first, in the top-level
% SuiteSparse folder.
%
% You may need to add the top-level lib folder (SPEX/lib, or SuiteSparse/lib
% if SPEX is inside SuiteSparse) to your LD_LIBRARY_PATH (DYLD_LIBRARY_PATH
% on the Mac).  See instructions in the spex_deps.m file.

% SPEX: (c) 2022, Chris Lourenco, Jinhao Chen,
% Lorena Mejia Domenzain, Timothy A. Davis, and Erick Moreno-Centeno.
% All Rights Reserved.
% SPDX-License-Identifier: GPL-2.0-or-later or LGPL-3.0-or-later


if (nargin < 1)
    run_demo = true ;
end

fprintf ('Compiling the SPEX mexFunctions for use:\n') ;

% Find all source files and add them to the src string
src = '';
path = './Source/';
files = dir('./Source/*.c');
[m n] = size(files);
for k = 1:m
    tmp = [' ', path, files(k).name];
    src = [src, tmp];
end
path = '../SPEX_Utilities/Source/';
files = dir('../SPEX_Utilities/Source/*.c');
[m n] = size(files);
for k = 1:m
    tmp = [' ', path, files(k).name];
    src = [src, tmp];
end

path = '../SPEX_LU/Source/';
files = dir('../SPEX_LU/Source/*.c');
[m n] = size(files);
for k = 1:m
    tmp = [' ', path, files(k).name];
    src = [src, tmp];
end

path = '../SPEX_Cholesky/Source/';
files = dir('../SPEX_Cholesky/Source/*.c');
[m n] = size(files);
for k = 1:m
    tmp = [' ', path, files(k).name];
    src = [src, tmp];
end
path = '../SPEX_Backslash/Source/';
files = dir('../SPEX_Backslash/Source/*.c');
[m n] = size(files);
for k = 1:m
    tmp = [' ', path, files(k).name];
    src = [src, tmp];
end

% Compiler flags
flags = 'CFLAGS=''-std=c99 -fPIC''';

% External libraries: GMP, MPRF, AMD, and COLAMD
[suitesparse_libdir, suitesparse_incdir, gmp_lib, gmp_include, mpfr_lib, mpfr_include] = spex_deps ;

% libraries:
if (isempty (suitesparse_libdir))
    suitesparse_libdir = ' ' ;
else
    suitesparse_libdir = [' -L' suitesparse_libdir ' '] ;
end
libs = [suitesparse_libdir ' -lamd -lcolamd -lsuitesparseconfig ' gmp_lib ' ' mpfr_lib ' -lm'] ;

% Path to headers
if (isempty (suitesparse_incdir))
    suitesparse_incdir = ' ' ;
else
    suitesparse_incdir = ['-I' suitesparse_incdir ' '] ;
end
if (isempty (gmp_include))
    gmp_include = ' ' ;
else
    gmp_include = [' -I' gmp_include ' '] ;
end
if (isempty (mpfr_include))
    mpfr_include = ' ' ;
else
    mpfr_include = [' -I' mpfr_include ' '] ;
end

includes = [ suitesparse_incdir ' -ISource/ -I../Include/ -I../../SuiteSparse_config -I../../COLAMD/Include -I../../AMD/Include -I../SPEX_Utilities/Source ' gmp_include  mpfr_include ] ;

% verbose = ' -v '
verbose = '' ;

% Generate the mex commands here
% having -R2018a here for function mxGetDoubles
m1 = ['mex ', verbose, ' -R2018a ', includes, ' spex_lu_mex_soln.c ' , src, ' ', flags, ' ', libs] ;
m2 = ['mex ', verbose, ' -R2018a ', includes, ' spex_cholesky_mex_soln.c ' , src, ' ', flags, ' ', libs];
m3 = ['mex ', verbose, ' -R2018a ', includes, ' spex_backslash_mex_soln.c ' , src, ' ', flags, ' ', libs];
m4 = ['mex ', verbose, ' -R2018a ', includes, ' spex_qr_mex_soln.c ' , src, ' ', flags, ' ', libs];

if (~isempty (verbose))
    fprintf ('%s\n', m1) ;
end

% Now, we evaluate each one
eval (m1) ;
eval (m2) ;
eval (m3) ;
eval (m4) ;

if (run_demo)
    % Test SPEX
    spex_mex_test ;
end

fprintf ('To use SPEX MATLAB Interface in future MATLAB sessions, add the following\n') ;
fprintf ('line to your startup.m file:\n') ;
fprintf ('   addpath (''%s'') ;\n', pwd) ;
fprintf ('Type ''doc startup'' for more info on how to use startup.m\n') ;
fprintf ('To run a demo, type:\n') ;
fprintf ('   echodemo spex_mex_demo ;\n') ;

