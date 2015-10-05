function [L,Plist,Nlist] = setup(N,domain,c)
%Assuming n = 2^k-1
l = 1;
Nl = N;
Nlist = Nl;
while Nl > 10
    l = l + 1;
    Nl = ((sqrt(Nl)-1)/2)^2;
    Nlist(l) = Nl;
end
L = l;
Nlist = fliplr(Nlist); % coarsest level: l = 1, finest level: l = L

Plist = {};
for l = L:-1:2
    Plist{l} = lin_interp(Nlist(l-1),Nlist(l));
    %Plist{l} = harmonic_interp(Nlist(l-1),Nlist(l),domain,c);
end

end
