function A = makematrix(domain,n,c)
% produce the matrix for div(c*grad(u))

x = linspace(domain{1}(1),domain{1}(2),n)';
y = linspace(domain{2}(1),domain{2}(2),n)';

h = (domain{1}(2)-domain{1}(1))/(n+1);

% x = domain(1,:)'; %column vector
% y = domain(2,:)'; %column vector

% n = size(x,1); % assume nx=ny
% h = x(2)-x(1); % assume uniform partitioning

Tj = @(x,yj) spdiags([-c(x+h/2,yj)...
    c(x+h/2,yj)+c(x-h/2,yj)+c(x,yj-h/2)+c(x,yj+h/2)...
    -c(x-h/2,yj)],[-1 0 1],n,n);  

cIj = @(x,yj) spdiags(c(x,yj+h/2),0,n,n);

Amd = zeros(n,n,n);
As = zeros(n,n,n-1);
for j = 1:n
    Amd(:,:,j) = Tj(x,y(j));
end

for j = 1:n-1
    As(:,:,j) = cIj(x,y(j));
end

A = blktridiag(Amd,-As,-As);

% T = blktridiag(2*(cx+cy),-cx,-cx,n);
% A = blktridiag(T,-cy*speye(n),-cy*speye(n),n);