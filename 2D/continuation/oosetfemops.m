function p=oosetfemops(p) % for 1D Neuman bc
fem=p.pdeo.fem; gr=p.pdeo.grid; 
x=getpte(p); x=x';
[p.mat.K,p.mat.M,~]=fem.assema(gr,x,1,1);

