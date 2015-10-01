clear all
n = 2^5-1; % (generalization needed)
h = 1/(n+1);
c = @(x,y) 5 + 2*cos(4*pi*x) + 2*cos(2*pi*y); 

x = h:h:1-h;
y = x;
domain = [x;y];

A = makematrix(domain,c);

% Make smoother
D = diag(diag(A));
L = tril(D-A);
Mrbgs = L*(D\L');
Mjac = D;

%% Defining system
g = -9.81;

b = h^2*g*ones(n^2,1);
% if n == 225
%     load('randb_p4.mat');
% else 
%     warning('n must be 225 if this b is to be used');
%     break
% end
tic
%% solving with multiple v-cycles
x = zeros(n^2,1);
m = 4;
resid_norm0 = sqrt(b'*b);
resid_norm = resid_norm0;
k = 0;
while resid_norm/resid_norm0 > 1e-4 && k<1000
    x = MV(A,b,Mjac,m,x);
    resid_norm = sqrt((A*x-b)'*(A*x-b));
    k = k + 1;
end
toc
k
%% Plotting
Uxy = reshape(x,n,n);
Ubar = zeros(n+2,n+2);
Ubar(2:end-1,2:end-1) = Uxy;
[X,Y] = meshgrid(linspace(0,1,n+2),linspace(0,1,n+2));
surf(X,Y,Ubar)
xlabel('x (c1)')
ylabel('y (c2)')
