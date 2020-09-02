function SPEX_backslash_install(run_demo)
% SPEX_backslash_INSTALL: install and test the MATLAB interface to SPEX_backslash.
% This function installs the SPEX backslash mexFunction for use by the m-file
% SPEX_backslash.m.
%
% Usage: SPEX_backslash_install
%
% Required Libraries: GMP, MPFR, AMD, COLAMD and SPEX.  If -lamd and -lcolamd are not
% available, install them with 'make install' first, in the top-level
% SuiteSparse folder.
%
% See also SPEX_Left_LU_backslash, SPEX_Cholesky_backslash

% SPEX: (c) 2019-2020, Chris Lourenco, Jinhao Chen, Erick Moreno-Centeno,
% Timothy A. Davis  All Rights Reserved.  See
% SPEX/License for the license.

if (nargin < 1)
    run_demo = true ;
end

fprintf ('Compiling the SPEX mexFunction for use in SPEX_backslash:\n') ;
fprintf ( 'This may take a few minutes, please wait\n');

% Install SPEX_Left_LU
cd ../SPEX_Left_LU/MATLAB/ ;
SPEX_Left_LU_install;

% Temporarily add SPEX_Left_LU to path
s = pwd;
addpath(s);

% Install SPEX_Chol
cd ../../SPEX_Cholesky/MATLAB ;
SPEX_Chol_install ;

% Temporarily add SPEX_Chol to path
s = pwd;
addpath(s);

% Return to this folder
cd ../../MATLAB ;


fprintf ('To use SPEX_backslash in future MATLAB sessions, add the following\n') ;
%TODO FIXME
fprintf ('line to your startup.m file:\n') ;
fprintf ('   addpath (''%s'') ;\n', pwd) ;
fprintf ('Type ''doc startup'' for more info on how to use startup.m\n') ;
fprintf ('To run a demo, type:\n') ;
fprintf ('   echodemo SPEX_backslash_demo ;\n') ;

