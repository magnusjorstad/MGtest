function showcoeff(c,domain)

x = linspace(domain{1}(1),domain{1}(2),50);
y = linspace(domain{2}(1),domain{2}(2),50);

[X,Y] = meshgrid(x,y);

surf(X,Y,c(X,Y))
end