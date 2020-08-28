%------------------------------------------------------------------------------
% IP-Chol/MATLAB/IP-Chol.m: Interface to IP Chol within MATLAB
%------------------------------------------------------------------------------

% IP_Chol: (c) 2020, Chris Lourenco, Erick Moreno-Centeno, Timothy A. Davis, 
% Texas A&M University.  All Rights Reserved.  See IP_Chol/License for the license.

%------------------------------------------------------------------------------

% TODO mimic SPEX backslash and Tim's changes
function varargout = IP_Chol(A,b,option)
% Purpose: Solve exactly the sparse linear system Ax = b where A and b are stored as
% doubles. A must be stored as a SPD sparse matrix. b must be stored as a dense
% set of right hand side vectors. b can be either 1 or multiple vector(s).
% In the matlab interface, the up-looking integer-preserving Chol is used
%
% ****WARNING****: If A is very large or dense, this function may crash.
%  
% USAGE:
% x = IP_Chol(A,b) returns the solution to Ax=b using default settings. A
% must be SPD
%
% x = IP_Chol(A,b,options) returns the solution to Ax=b with user defined
% settings. The options settings can be obtained from option =
% IP_get_options then changed from there
%

if exist('option') == 0          % Did the user pass in options?
    option = IP_get_options;   % Set defaults
end

% Check for int overflow. If the max value of A or b exceeds this value
% the internal routines can not expect int input.
if (max(max(abs(A))) > 2000000000000) 
    option.int = 0;
end
if (max(abs(A)) > 2000000000000)
    option.intb = 0;
end

% Check if the input matrix is stored as sparse. If not, SPEX LU expects
% sparse input, so convert to sparse.
if (issparse(A) == 0)
    A = sparse(A);
end

% If the user indicates that the input is integral, check if it is actually integral
if (option.int > 0)
    A2 = floor(A);
    if (normest(A2-A) > 1e-12)
        option.int = 0;
    end
    clear A2;
end

if (option.intb > 0)
    b2 = floor(b);
    if (normest(b2-b) > 1e-12)
        option.intb = 0;
    end
    clear b2;
end

% Matrix has cleared input checks. Now test to see if cholesky
% factorization runs to completion.

try chol(A);
% MATLAB flags the matrix as SPD, attempt to use IP_Chol

if (nargout == 1) % x = A\b
    varargout{1} = IP_mex_soln(A,b,option);
else
fprintf('Incorrect number of output arguments. Please type help IP_Chol\n')
end

catch
    fprintf('Input matrix is not SPD. To solve your system, please use SPEX LU\n')
end

end
