clear all
n = 2^6-1; % (generalization needed)
h = 1/(n+1);
c = @(x,y) 1 + 100*(heaviside(x-0.2)-heaviside(x-0.8));

domain{1} = [0 1];
domain{2} = [0 1];
%showcoeff(c,domain);

A = makematrix(domain,n,c);

%% RHS
g = -9.81;
b = h^2*g*ones(n^2,1);
% if n == 225
%     load('randb_p4.mat');
% else 
%     warning('n must be 2^4-1=225 if this b is to be used');
%     break
% end


%% solving with multigrid
[x,k] = MG(A,b,domain,c);


%% Plotting
Uxy = reshape(x,n,n);
Ubar = zeros(n+2,n+2);
Ubar(2:end-1,2:end-1) = Uxy;
[X,Y] = meshgrid(linspace(0,1,n+2),linspace(0,1,n+2));
surf(X,Y,Ubar)
xlabel('x (c1)')
ylabel('y (c2)')
