function [x,k] = MG(A,b,domain,c)

N = size(A,1);
m = 4;

% Make smoother
D = diag(diag(A));
[L,U] = lu(A);
L = L - D;
Mgs = tril(A);
Mjac = D;

% Set up prolongation operators for all grids
tic
[L,Plist,Nlist] = setup(N,domain,c);
fprintf(1,'Setup time: ')
toc


tic
resid_norm0 = sqrt(b'*b);
resid_norm = resid_norm0;
x = zeros(Nlist(L),1);
k = 0;
while resid_norm/resid_norm0 > 1e-4 && k<1000
    x = MV(A,b,Mgs,m,x,Plist,Nlist,L,L);
    resid_norm = sqrt((A*x-b)'*(A*x-b));
    k = k + 1;
end
fprintf(1,'MG solve time: ')
toc
end