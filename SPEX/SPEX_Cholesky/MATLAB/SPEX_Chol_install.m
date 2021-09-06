function SPEX_Chol_install(run_demo)
%SPEX_Chol_INSTALL: install and test the MATLAB interface to SPEX_backslash.
% This function installs the SPEX Chol mexFunction for use by the m-file
% SPEX_Chol_backslash.m.
%
% Usage: SPEX_Chol_install
%
% Required Libraries: GMP, MPFR, AMD, COLAMD, SPEX.  If -lamd and -lcolamd are not
% available, install them with 'make install' first, in the top-level
% SuiteSparse folder.
%
% See also SPEX_Chol_backslash, SPEX_Chol_test, SPEX_Chol_demo.

% SPEX_Cholesky: (c) 2020, Chris Lourenco, United States Naval Academy, 
% Erick Moreno-Centeno, Timothy A. Davis, Jinhao Chen, Texas A&M University.  
% All Rights Reserved.  See SPEX_Cholesky/License for the license.
n
if (nargin < 1)
    run_demo = true ;
end

fprintf ('Compiling the SPEX Chol mexFunction for use in SPEX_Chol_backslash:\n') ;

% Find all source files and add them to the src string
src = '';
path = './Source/';
files = dir('./Source/*.c');
[m n] = size(files);
for k = 1:m
    tmp = [' ', path, files(k).name];
    src = [src, tmp];
end
path = '../Source/';
files = dir('../Source/*.c');
[m n] = size(files);
for k = 1:m
    tmp = [' ', path, files(k).name];
    src = [src, tmp];
end
path = '../../SPEX_Util/Source/';
files = dir('../../SPEX_Util/Source/*.c');
[m n] = size(files);
for k = 1:m
    tmp = [' ', path, files(k).name];
    src = [src, tmp];
end

% Compiler flags
flags = 'CFLAGS=''-std=c99 -fPIC''';

% External libraries: GMP, MPRF, AMD, and COLAMD
libs = '-L../../lib -lgmp -lmpfr -lamd -lcolamd -lsuitesparseconfig' ;

% Path to headers
includes = '-ISource/ -I../Source/ -I../Include/ -I../../SPEX_Util/Include -I../../../SuiteSparse_config -I../../../COLAMD/Include -I../../../AMD/Include -I../../SPEX_Util/Source';

% verbose = ' -v '
verbose = '' ;

% Generate the mex commands here
% having -R2018a here for function mxGetDoubles
m1 = ['mex ', verbose, ' -R2018a ', includes, ' SPEX_Chol_mex_soln.c ' , src, ' ', flags, ' ', libs];

if (~isempty (verbose))
    fprintf ('%s\n', m1) ;
end

% Now, we evaluate each one
eval (m1) ;

if (run_demo)
    % Test SPEX_backslash.
    SPEX_Chol_test ;
end

fprintf ('To use SPEX_Chol_backslash in future MATLAB sessions, add the following\n') ;
fprintf ('line to your startup.m file:\n') ;
fprintf ('   addpath (''%s'') ;\n', pwd) ;
fprintf ('Type ''doc startup'' for more info on how to use startup.m\n') ;
fprintf ('To run a demo, type:\n') ;
fprintf ('   echodemo SPEX_Chol_demo ;\n') ;

