function v = MV(A,f,M,m,v0,Rf,IR)
% Recursive multigrid V-cycle. 2D grid, n blocks of size n, n=2^k-1
nf = sqrt(length(v0));
v = v0;

if nf < 10
    v = A\f; %Exact solution
else
    
    %Pre-smoothing m times
    for i = 1:m
        v  = M\((M-A)*v+f);       
    end
    
    % Prolongation operator. Restriction is defined as R = P'/4
    P = prolongation(Nc,Nf);
    
    
    r = A*v-f;
    Ac = 4\P'*A*P;
    rc = 4\P'*r;
    Mc = 4\P'*M*P;
    
    vhat = MV(Ac,rc,Mc,m,zeros(size(Ac,1),1));
    
    nc = size(vhat,1);
    P = prolongation(nc,nf);
    v = v + P*vhat;
    
%     
%     %Mappings, Rf:fine->coarse, I:coarse->fine
%     nc = 2^(log2(nf+1)-1)-1;
%     ncc = max(2^(log2(nc+1)-1)-1,1);
%     Il = sparse([eye(ncc) zeros(ncc,nc-ncc)]);
%     IL = kron(Il,Il);
%     
%     If = 4*Rf';
%     %Rf*(f-A*v)
%     %Compute error on coarser grid
%     vhat = MV(Rf*A*If,Rf*(f-A*v),Rf*M*If,m,zeros(size(Rf*f)),IL*Rf*IR,IL');
%     
%     %Coarse grid correction
%     v = v + If*vhat;
%     

    %Post-smoothing m times
    for i = 1:m
        v  = M\((M-A)*v+f);
    end
end
    

    