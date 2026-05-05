function p=oosetfemops(p) % for 1D Neuman bc
fem=p.pdeo.fem; gr=p.pdeo.grid; 
[p.mat.K,p.mat.M,~]=p.pdeo.fem.assema(gr,1,1,1); 
x=getpte(p); x=x';
p.mat.Kx=convection(fem,gr,1/x); % first derivative respect to x (needed for the laplacian in polar coordinates)