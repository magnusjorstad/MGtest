function v = Mcycle(A,f,M,m,v0,Plist,Nlist,l,L)
% Recursive multigrid V-cycle (or w-cycle). 
% 2D grid, n blocks of size n, n=2^k-1

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
    
    % Restriction
    [rc,Ac,Mc] = res(f-A*v,A,M,P);

    %Coarse grid solution
    vhat = Mcycle(Ac,rc,Mc,m,zeros(Nc,1),Plist,Nlist,l-1,L);
    
    if 0 % w-cycle
        vhat = Mcycle(Ac,rc,Mc,m,vhat,Plist,Nlist,l-1,L);
    end
    
    % Prolongation
    v = v + P*vhat;

    %Post-smoothing m times
    for i = 1:m
        v  = M\((M-A)*v+f);
    end
end
    

    