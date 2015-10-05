function [rc,Ac,Mc] = res(r,A,M,P)
    Ac = 4\P'*A*P;
    rc = 4\P'*r;
    Mc = 4\P'*M*P;
end