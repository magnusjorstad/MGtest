function v = MV(A,f,M,m,v0,Plist,Nlist,l,L)
% Recursive multigrid V-cycle. 2D grid, n blocks of size n, n=2^k-1

if l == 1
    v = A\f; %Exact solution
else
    v = v0;
    Nc = Nlist(l-1);
    %Pre-smoothing m times
    for i = 1:m
        v  = M\((M-A)*v+f);       
    end
    
    % Prolongation operator. P: G(l-1) -> G(l), R := P'/4
    P = Plist{l};
    
    % Coarse grid parameters and variables
    r = f-A*v;
    Ac = 4\P'*A*P;
    rc = 4\P'*r;
    Mc = 4\P'*M*P;
    
    vhat = MV(Ac,rc,Mc,m,zeros(Nc,1),Plist,Nlist,l-1,L);
    
    v = v + P*vhat;

    %Post-smoothing m times
    for i = 1:m
        v  = M\((M-A)*v+f);
    end
end
    

    