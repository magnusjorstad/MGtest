clear all
n = 3;
c1 = 10; c2 = 1; % anisotropic
T = blktridiag(2*(c1+c2),-c1,-c1,n);
A = blktridiag(T,-c2*speye(n),-c2*speye(n),n);
D = diag(diag(A));
L = tril(D-A);
Mrbgs = L*D\L';
Mjac = D;

%% Defining system
g = -9.81;
h = 1/(n+1);
b = h^2*g*ones(n^2,1);

%U = A\b;

%% solving with multiple v-cycles
U = zeros(n^2,1);
m = 2;
resid_norm0 = sqrt(b'*b);
resid_norm = resid_norm0;
k = 0;
while resid_norm/resid_norm0 > 1e-4 && k<100
    x = MV(A,b,Mjac,m,U);
    resid_norm = sqrt((A*x-b)'*(A*x-b));
    k = k + 1;
end

%% Plotting
Uxy = reshape(U,n,n);
Ubar = zeros(n+2,n+2);
Ubar(2:end-1,2:end-1) = Uxy;
[X,Y] = meshgrid(linspace(0,1,n+2),linspace(0,1,n+2));
surf(X,Y,Ubar)
xlabel('x (c1)')
ylabel('y (c2)')