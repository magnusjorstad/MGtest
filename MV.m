function v = MV(A,f,M,m,v0)
% Recursive multigrid V-cycle. 2D grid, n blocks of size n, n=2^k-1
Nf = length(v0);
Nc = ((sqrt(Nf)-1)/2)^2;
v = v0;

if Nf < 10
    v = A\f; %Exact solution
else
    
    %Pre-smoothing m times
    for i = 1:m
        v  = M\((M-A)*v+f);       
    end
    
    % Prolongation operator. Restriction is defined as R = P'/4
    P = prolongation(Nc,Nf);
    
    % Coarse grid parameters and variables
    r = A*v-f;
    Ac = 4\P'*A*P;
    rc = 4\P'*r;
    Mc = 4\P'*M*P;
    
    vhat = MV(Ac,rc,Mc,m,zeros(Nc,1));
    
    %nc = size(vhat,1);
    %P = prolongation(nc,nf);
    v = v + P*vhat;

    %Post-smoothing m times
    for i = 1:m
        v  = M\((M-A)*v+f);
    end
end
    

    