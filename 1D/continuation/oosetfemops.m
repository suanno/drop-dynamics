function p=oosetfemops(p) % for 1Dpbc, with x-dep. K 
fem=p.pdeo.fem; gr=p.pdeo.grid; 
[K,M,~]=p.pdeo.fem.assema(p.pdeo.grid,1,1,1); 
Kx=convection(fem,gr,1); p.mat.Kx=filltrafo(p,Kx); % Convection matrix needed for phase condition (translational invariance)
p.mat.M0=p.mat.fill'*M; % we need M0 to transform the nonlinearity 
p.mat.K=filltrafo(p,K); p.mat.M=filltrafo(p,M); % transform of K and M 