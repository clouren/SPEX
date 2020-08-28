% This function installs the integer-preserving Cholesky matlab routines. 
% It allows the use of m files IP_Chol.m
% Please run this command prior to attempting to use any IP Chol routines
% Usage: IP_install
% Required Libraries: GMP, MPFR, AMD, COLAMD, SPEX LU

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

path = '../SPEX_LU-master/SPEX_LU/Source/';
files = dir('../SPEX_LU-master/SPEX_LU/Source/*.c');
[m n] = size(files);
for k = 1:m
    tmp = [' ', path, files(k).name];
    src = [src, tmp];
end


% Compiler flags
flags = 'CFLAGS=''-std=c99 -fPIC''';

% External libraries
libs = '-lgmp -lmpfr -lamd -lcolamd -lsuitesparseconfig';

% Path to headers
includes = '-ISource/ -I../Source/ -I../Include/ -I../';

includes = ['-I../Include -I../Source -ISource -I../SPEX_LU-master/SuiteSparse_config ' ...
'-I../SPEX_LU-master/COLAMD/Include -I../SPEX_LU-master/AMD/Include -I../SPEX_LU-master/SPEX_LU/Include ' ...
'-I../SPEX_LU-master/SPEX_LU/Source -I../SPEX_LU-master/SPEX_LU/Lib'];

% Generate the mex commands here
% having -R2018a here for function mxGetDoubles
m1 = ['mex -R2018a ', includes, ' IP_mex_soln.c ' , src, ' ', flags, ' ', libs];

% Now, we evaluate each one
eval(m1);

fprintf('\nMex files installed, now we test\n')

% Efficient testing
IP_test;
