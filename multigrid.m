function [error,k] = multigrid(n,m,omega)
T = blktridiag(4,-1,-1,n);
A = blktridiag(T,-speye(n),-speye(n),n);
D = diag(diag(A));
L = tril(D-A);

if nargin < 3
omega = 1;
end
Mssor = 0.25*(D-omega*L)*(D-omega*L')/(omega*(2-omega));


%% Defining system
x_sol= rand(n^2,1); 
load x_sol
if n == 15, x_sol = x15;end
if n == 31, x_sol = x31;end
if n == 63, x_sol = x63;end
if n == 127, x_sol = x127;end
if n == 255, x_sol = x255;end

b = A*x_sol;
x = zeros(n^2,1);
error0 = sqrt(x_sol'*A*x_sol);



%% solving with multiple v-cycles
error = error0;
k = 0;
Rf = sparse(restriction(n));
nc = 2^(log2(n+1)-1)-1;
Ir = sparse([eye(nc);zeros(n-nc,nc)]);
IR = kron(Ir,Ir);
while error/error0 > 1e-6 && k<100
    x = MV(A,b,Mssor,m,x,Rf,IR);
    error = sqrt((x-x_sol)'*A*(x-x_sol));
    k = k + 1;
end

end