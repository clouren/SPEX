function SPEX_mex_install(run_demo)
% SPEX_mex_INSTALL: install and test the MATLAB interface to SPEX MATLAB functions.
%
% Usage: SPEX_mex_install
%
% Required Libraries: GMP, MPFR, AMD, COLAMD, SPEX.  If -lamd and -lcolamd are not
% available, install them with 'make install' first, in the top-level
% SuiteSparse folder.

% SPEX: (c) 2022, Chris Lourenco, United States Naval Academy, 
% Jinhao Chen, Lorena Mejia Domenzain, Erick Moreno-Centeno and Timothy A. Davis, 
% Texas A&M University.  All Rights Reserved.  
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
libs = '-L/home/grads/l/lorena.m.d/Documents/RESEARCH/SPEX/lib -lamd -lcolamd -lsuitesparseconfig -L/home/grads/l/lorena.m.d/spack/opt/spack/linux-ubuntu18.04-x86_64_v4/gcc-7.5.0/gmp-6.2.1-kbb7qvlehyep4sbupog6a4ugtrcjlte4/lib -lgmp -lm -L/home/grads/l/lorena.m.d/spack/opt/spack/linux-ubuntu18.04-x86_64_v4/gcc-7.5.0/mpfr-4.1.0-qdbpzu652zw5zkmmngvhrjxg6iokbsls/lib -lmpfr -lm' ;

% Path to headers
includes = '-ISource/ -I../SPEX_Backslash/Include/ -I../../SuiteSparse_config -I../../COLAMD/Include -I../../AMD/Include -I../SPEX_Utilities/Source -I/home/grads/l/lorena.m.d/spack/opt/spack/linux-ubuntu18.04-x86_64_v4/gcc-7.5.0/gmp-6.2.1-kbb7qvlehyep4sbupog6a4ugtrcjlte4/include -I/home/grads/l/lorena.m.d/spack/opt/spack/linux-ubuntu18.04-x86_64_v4/gcc-7.5.0/mpfr-4.1.0-qdbpzu652zw5zkmmngvhrjxg6iokbsls/include';

% verbose = ' -v '
verbose = '' ;

% Generate the mex commands here
% having -R2018a here for function mxGetDoubles
m1 = ['mex ', verbose, ' -R2018a ', includes, ' spex_lu_mex_soln.c ' , src, ' ', flags, ' ', libs];
m2 = ['mex ', verbose, ' -R2018a ', includes, ' spex_cholesky_mex_soln.c ' , src, ' ', flags, ' ', libs];
m3 = ['mex ', verbose, ' -R2018a ', includes, ' spex_backslash_mex_soln.c ' , src, ' ', flags, ' ', libs];

if (~isempty (verbose))
    fprintf ('%s\n', m1) ;
end

% Now, we evaluate each one
eval (m1) ;
eval (m2) ;
eval (m3) ;

if (run_demo)
    % Test SPEX
    SPEX_mex_test ;
end

fprintf ('To use SPEX MATLAB Interface in future MATLAB sessions, add the following\n') ;
fprintf ('line to your startup.m file:\n') ;
fprintf ('   addpath (''%s'') ;\n', pwd) ;
fprintf ('Type ''doc startup'' for more info on how to use startup.m\n') ;
fprintf ('To run a demo, type:\n') ;
fprintf ('   echodemo SPEX_mex_demo ;\n') ;
