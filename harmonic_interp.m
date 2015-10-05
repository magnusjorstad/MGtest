function P = harmonic_interp(Nc,Nf,domain,c)
% P maps from coarse to fine
% Assumes uniform square grid with nx=ny

nc = sqrt(Nc);
nf = sqrt(Nf);

x = linspace(domain{1}(1),domain{1}(2),nc);
y = linspace(domain{2}(1),domain{2}(2),nc);
hf = (domain{1}(2)-domain{1}(1))/(nf+1); % hf = (b-a)/(nf+1)
P = zeros(Nf,Nc);

for iy = 1:nc
    for ix = 1:nc
    ic = (2*iy-1)*nf + 2*ix;
    
    % %%%
    % %0%
    % %%%
    P(ic,nc*(iy-1)+ix) = 1;
    
    % %%%
    % 0%%
    % %%%
    P(ic-1,nc*(iy-1)+ix) = c(x(ix)-hf/2,y(iy)) / ...
        (c(x(ix)-hf/2,y(iy)) + c(x(ix)-2*hf/3,y(iy)));
    
    % %%%
    % %%0
    % %%%
    P(ic+1,nc*(iy-1)+ix) = c(x(ix)+hf/2,y(iy)) / ...
        (c(x(ix)+hf/2,y(iy)) + c(x(ix)+2*hf/3,y(iy)));
        
    % %%%
    % %%%
    % %0%
    P(ic-nf,nc*(iy-1)+ix) = c(x(ix),y(iy)-hf/2) / ...
        (c(x(ix),y(iy)-hf/2) + c(x(ix),y(iy)-2*hf/3));
    
    % %0%
    % %%%
    % %%%
    P(ic+nf,nc*(iy-1)+ix) = c(x(ix),y(iy)+hf/2) / ...
        (c(x(ix),y(iy)+hf/2) + c(x(ix),y(iy)+2*hf/3));
    
    % %%%
    % %%%
    % 0%%
    P(ic-nf-1,nc*(iy-1)+ix) = ...
        c(x(ix)-hf/2,y(iy)-hf)*P(ic-nf,nc*(iy-1)+ix) / ...
        (c(x(ix)-hf/2,y(iy)-hf) + c(x(ix)-3*hf/2,y(iy)-hf)) + ...
        c(x(ix)-hf,y(iy)-hf/2)*P(ic-1,nc*(iy-1)+ix) / ...
        (c(x(ix)-hf,y(iy)-hf/2) + c(x(ix)-hf,y(iy)-3*hf/2));
    
    % %%%
    % %%%
    % %%0
    P(ic-nf+1,nc*(iy-1)+ix) = ...
        c(x(ix)+hf/2,y(iy)-hf)*P(ic-nf,nc*(iy-1)+ix) / ...
        (c(x(ix)+hf/2,y(iy)-hf) + c(x(ix)+3*hf/2,y(iy)-hf)) + ...
        c(x(ix)+hf,y(iy)-hf/2)*P(ic+1,nc*(iy-1)+ix) / ...
        (c(x(ix)+hf,y(iy)-hf/2) + c(x(ix)+hf,y(iy)-3*hf/2));
    
    % 0%%
    % %%%
    % %%%
    P(ic-nf+1,nc*(iy-1)+ix) = ...
        c(x(ix)-hf/2,y(iy)+hf)*P(ic+nf,nc*(iy-1)+ix) / ...
        (c(x(ix)-hf/2,y(iy)+hf) + c(x(ix)-3*hf/2,y(iy)+hf)) + ...
        c(x(ix)-hf,y(iy)+hf/2)*P(ic-1,nc*(iy-1)+ix) / ...
        (c(x(ix)-hf,y(iy)+hf/2) + c(x(ix)-hf,y(iy)+3*hf/2));
    
    
    % %%0
    % %%%
    % %%%
    P(ic-nf-1,nc*(iy-1)+ix) = ...
        c(x(ix)+hf/2,y(iy)+hf)*P(ic+nf,nc*(iy-1)+ix) / ...
        (c(x(ix)+hf/2,y(iy)+hf) + c(x(ix)+3*hf/2,y(iy)+hf)) + ...
        c(x(ix)+hf,y(iy)+hf/2)*P(ic+1,nc*(iy-1)+ix) / ...
        (c(x(ix)+hf,y(iy)+hf/2) + c(x(ix)+hf,y(iy)+3*hf/2)); 
    
    
    end
end
P = sparse(P);
