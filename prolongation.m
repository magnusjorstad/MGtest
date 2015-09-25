function P = prolongation(Nc,Nf)
% Map from coarse to fine
nc = sqrt(Nc);
nf = sqrt(Nf);

Pbar = zeros(nf,nc);

for i = 1:nc
    Pbar((2*i-1):(2*i+1),i) = [.5;1;.5];
end


P = kron(Pbar,Pbar);
