clear all
n = 19;
T = blktridiag(4,-1,-1,n);
A = blktridiag(T,-speye(n),-speye(n),n);
D = diag(diag(A));
L = tril(D-A);

%% Defining system
% % Random
% x_sol= rand(n^2,1); 
% b = A*x_sol;
% x_0 = zeros(n^2,1);

% Gravitation
h = 1/(n+1);
g = -9.81;
b = h^2*g*ones(n^2,1);
midpoint = n*floor(n/2)+ceil(n/2);
U = A\b;
U(midpoint)
Uxy = reshape(U,n,n);
Ubar = zeros(n+2,n+2);
Ubar(2:end-1,2:end-1) = Uxy;
[X,Y] = meshgrid(linspace(0,1,n+2),linspace(0,1,n+2));
surf(X,Y,Ubar)
break
%% Preconditioners
I = speye(n^2); % gives pure CG

Mjacobi = 4*I;

omega = 3/2;
Mssor = (D-omega*L)*(D-omega*L')/(omega*(2-omega))/4;
    
%% Solving system
tic;
[x, k, gamma,error] = pCG(x_0, A, Mjacobi, b,x_sol, 20000, 1e-9, @(Z,o) Z*o, @(Z,o) Z\o);
toc;
gm = mean(gamma);
genK = ((1+gm/2)/(1-gm/2))^2;
meanrate = mean(gamma(ceil(k/5):end));

k
% figure(1);hold off;
% plot(log(error),'b');
% hold on
% plot(1:k,-(1:k),'r');
