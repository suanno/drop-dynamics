function out = bra(p, u)
H = u(1:p.nu); % Exclude the parameters of the pde
par = u(p.nu+1:end); % Parameters
R=getpte(p); R=R';
e1 = ones(1,length(R))';

c0 = sum(p.mat.M*(H))/p.vol;
Hmax = max(H);
Hout = 1+(wetting_potential(Hmax)-wetting_potential(1))/((Hmax-1)*second_deriv_wetting_potential(1));
Hout_ = H(end);
%Hout = Hout_; % The estimate of Hout from Hmax is wrong! Problem as it
%makes the result depend on the box size as integrals of Hhat do not converge with
%box size going to infinity!
Hhat = H-Hout_*e1;

% %dhomegaout = par(3);
% %TODO: In 2D there is a factor R in the integrand, TO FIX!
% omegahat = wetting_potential(H)-wetting_potential(Hout)*e1;
% Omegahat = sum(p.mat.M*(omegahat));
% m = sum(p.mat.M*Hhat);
% I = sum(p.mat.M*(Hout^3.*Hhat./H.^3));
% K = sum(p.mat.M*(Hhat.^2./H.^3));

% Psi^nu
% - Solve ODE for Psi1
%G = (3*Hout_./H-2);
c = (H.^3).^(-1);
a = 1 ./ ((1e-12+R.^2) .* (H.^3));    % 1e-12 To eliminate the divergence at r=0
%Hr = diff(H)./diff(R);
%Hr = [0; Hr];
%frhs = -G.*(Hr./(H.^3));
frhs = -diff((H-Hout_)./H.^3)./diff(R);
frhs = [0; frhs];
fem=p.pdeo.fem;
gr=p.pdeo.grid;
[Kpsi,Mpsi,Fpsi] = fem.assema(gr,c,a,frhs);
% Robin bc at x=L 
A = Kpsi + Mpsi;
phi = (1 + sqrt(5))/2;
alpha_ = c(end)*(1 - phi)/R(end);
A(end,end) = A(end,end) - alpha_;
Psi1 = A\Fpsi;
Psi1lim = Psi1(end);

% - Solve ODE for Psi2
%G = (Hout_.^3)./H;
%frhs = -G.*(Hr./(H.^3));
frhs = (1/3)*diff(Hout_^3./H.^3)./diff(R);
frhs = [0; frhs];
[Kpsi,Mpsi,Fpsi] = fem.assema(gr,c,a,frhs);
% Robin bc at x=L 
A = Kpsi + Mpsi;
A(end,end) = A(end,end) - alpha_;  % see sign convention below
Psi2 = A\Fpsi;
%Psi2lim = Psi2(end);


% Drop radius
x0 = R(1);  H0 = H(1);
x1 = R(2);  H1 = H(2);
x2 = R(3);  H2 = H(3);
% Parabolic interpolant through the first three points
fit_coeff = polyfit([x0 x1 x2], [H0 H1 H2], 2);
% Coefficient of the parabolic approximation h(r) = Hmax + c*r^2
d2h0 = 2*fit_coeff(1);
c = 0.5*d2h0;      % = p(1)
% Compute where the parabola intersects h=hout
rcrit = sqrt((Hout_ - Hmax)/c);

%Contribution of B0in
% - v1 term
K_noB = 2*pi*sum(p.mat.M*(((Hhat.^2)./(H.^3)).*R));
RPsi1R = diff(R.*Psi1)./diff(R);
RPsi1R = [0; RPsi1R];
K_Bcorr = 3*pi*sum(p.mat.M*(((Hhat)./(H.^3)).*RPsi1R));
% - gradW term
I_noB = 2*pi*sum(p.mat.M*(((Hhat.*Hout_.^3)./(H.^3)).*R));
RPsi2R = diff(R.*Psi2)./diff(R);
RPsi2R = [0; RPsi2R];
I_Bcorr = 3*pi*sum(p.mat.M*(((Hhat)./(H.^3)).*RPsi2R));

%Integrals
Omega=2*pi*sum(p.mat.M*(R.*(wetting_potential(H)-wetting_potential(Hout_))));
V = 2*pi*sum(p.mat.M*(R.*Hhat));
I=2*pi*sum(p.mat.M*(((Hhat)./(H.^3)).*(R.*Hout_.^3-3*0.5*RPsi2R)));
K=2*pi*sum(p.mat.M*((Hhat./(H.^3)).*(Hhat.*R+3*0.5*RPsi1R)));
%dhomega = deriv_wetting_potential(H);
%dhOmega = sum(p.mat.M0*(dhomega));

out = [Omega;K;I;V; Hmax; Hout; K_noB; K_Bcorr; I_noB; I_Bcorr; Hout_; Psi1lim; c0; rcrit];