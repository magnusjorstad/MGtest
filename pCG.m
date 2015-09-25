function [x, k, gamma,error] = pCG(x0, A, M, b,xsol, mit, stol, bbA, bbM)
% Synopsis:
% x0: initial point
% A: Matrix A of the system Ax=b
% M: Preconditioning Matrix
% mit: Maximum number of iterations
% stol: energy norm tolerance
% bbA: Black Box that computes the matrix-vector product for A * u
% bbC: Black Box that computes: z = M\r
% x: Estimated solution point
% k: Number of iterations done 
%
% Example:
% tic;[x, t] = cgp(x0, S, speye(1), b, 3000, 10^-8, @(Z, v) Z*v, @(Z, v) Z\v);toc
% Elapsed time is 0.550190 seconds.

if ( nargin < 8 ), error('Not enough input arguments. Try help.'); end;
if ( isempty(A) ), error('Input matrix A must not be empty.'); end;
if ( isempty(M) ), error('Input preconditioner matrix M must not be empty.'); end;

x = x0;
p = 0;
k = 0;
e0 = x-xsol;e = e0;gamma = 0;error(1) = e0'*A*e0;


r = b - bbA(A, x0); % <--- r = b - A * x0;
while ( (e'*A*e)/(e0'*A*e0) > stol ),
    z = bbM(M, r); % <--- z = M \ r;
    k = k + 1;
    if ( k == mit ), warning('GCP:MAXIT', 'mit reached, no conversion.'); return; end;
    
    if ( k == 1 )
        p = z;
        rho = r'*z;
    else
        rho1 = rho;
        rho = r'*z;
        beta = rho/rho1;
        p = z + beta * p;
    end;
    Ap = bbA(A, p); % <--- Ap = A * p;
    a = rho / (p'*Ap);
    x = x + a * p;
    r = r - a * Ap;
    
    eprev = e;
    e = x-xsol;error(k+1) = e'*A*e;
    gamma(k) = sqrt((e'*A*e)/(eprev'*A*eprev));
end;