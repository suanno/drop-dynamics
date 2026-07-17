function out = bra(p, u)
H = u(1:p.nu); % Exclude the parameters of the pde
par = u(p.nu+1:end); % Parameters
R=getpte(p); R=R';
e1 = ones(1,length(R))';

c0 = sum(p.mat.M*(H))/p.vol;
Hmax = max(H);
Hout = 1+(wetting_potential(Hmax)-wetting_potential(1))/((Hmax-1)*second_deriv_wetting_potential(1));
Hhat = H-Hout*e1;

% %dhomegaout = par(3);
% %TODO: In 2D there is a factor R in the integrand, TO FIX!
% omegahat = wetting_potential(H)-wetting_potential(Hout)*e1;
% Omegahat = sum(p.mat.M*(omegahat));
% m = sum(p.mat.M*Hhat);
% I = sum(p.mat.M*(Hout^3.*Hhat./H.^3));
% K = sum(p.mat.M*(Hhat.^2./H.^3));

% - Solve ODE for Psi1
G = (3*Hout./H-2);
c = (H.^3).^(-1);
a = 1 ./ ((1e-12+R.^2) .* (H.^3));    % 1e-12 To eliminate the divergence at r=0
Hr = diff(H)./diff(R);
Hr = [0; Hr];
frhs = -G.*(Hr./(H.^3));

fem=p.pdeo.fem;
gr=p.pdeo.grid;
[Kpsi,Mpsi,Fpsi] = fem.assema(gr,c,a,frhs);
Psi1 = (Kpsi + Mpsi)\Fpsi;
% - Solve ODE for Psi2
G = (Hout.^3)./H;
frhs = -G.*(Hr./(H.^3));

[Kpsi,Mpsi,Fpsi] = fem.assema(gr,c,a,frhs);
Psi2 = (Kpsi + Mpsi)\Fpsi;


%Contribution of B0in
% - v1 term
I1v1 = 2*sum(p.mat.M*(((Hhat.^2)./(H.^3)).*R));
RPsi1R = diff(R.*Psi1)./diff(R.*Psi1);
RPsi1R = [0; RPsi1R];
I2v1 = 2*sum(p.mat.M*(((Hhat)./(H.^3)).*RPsi1R));
% - gradW term
I1gradW = (2/3)*sum(p.mat.M*(((Hhat.*Hout.^3)./(H.^3)).*R));
RPsi2R = diff(R.*Psi2)./diff(R.*Psi2);
RPsi2R = [0; RPsi2R];
I2gradW = 2*sum(p.mat.M*(((Hhat)./(H.^3)).*RPsi2R));

out = [Hmax; Hout; I1v1; I2v1; I1gradW; I2gradW];