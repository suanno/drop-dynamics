%% demo for pBC 1D, clear workspace 
close all; keep pphome; 
%% cell 1: init
p=[]; par=[0.5 0 0]; % concentration, chemical potential (lagrange mult mass)
p=chinit(p,12.5,2000,par); p.sw.qjac=0; %p.sw.verb=2;
%% Continuation parameters
p.nc.ilam = [2, 3];
p.nc.nq=1;
p.nc.lammax=20; p.sol.ds=0.1; p.nc.dsmax=0.5;
%% First branch continuation
p=setfn(p,'tr'); p=findbif(p,30);
%% Switch to droplet branch
p.sol.ds=+0.01;
p=swibra('tr','bpt1','b1',+0.01); p=cont(p,30); 
p.sw.bifcheck=0;
p.sw.foldcheck=1;

%%
H = p.u(1:p.nu); % Exclude the parameters of the pde
Hout = min(H);
Hmax = max(H);
rval=getpte(p); rval=rval';

% Compute drop radius
x0 = rval(1);  H0 = H(1);
x1 = rval(2);  H1 = H(2);
x2 = rval(3);  H2 = H(3);
% Parabolic interpolant through the first three points
fit_coeff = polyfit([x0 x1 x2], [H0 H1 H2], 2);
% Coefficient of the parabolic approximation h(r) = Hmax + c*r^2
d2h0 = 2*fit_coeff(1);
c = 0.5*d2h0;      % = p(1)
% Compute where the parabola intersects h=hout
rcrit = sqrt((Hout - Hmax)/c);


% Solve ode for Psi 1
g = H - Hout;
eta=1;
Qin = H.^3/(3*eta);

c = Qin.^(-1);
a = 1 ./ ((1e-12+rval.^2) .* Qin);    % To eliminate the divergence

fac  = 1 - 3*g./H;
Hr = diff(H)./diff(rval);
Hr = [0; Hr];
frhs = -(fac .* Hr) ./ Qin;

fem=p.pdeo.fem;
gr=p.pdeo.grid;
[Kpsi,Mpsi,Fpsi] = fem.assema(gr,c,a,frhs);
psi1 = (Kpsi + Mpsi)\Fpsi;

% Psi 2

Qin  = H.^3/(3*eta);
Qout = Hout.^3/(3*eta);
% c(r) and a(r) are the same! Only frhs is different

frhs = -3 .* (Hout.^3 ./ H.^4) .* Hr;

[Kpsi,Mpsi,Fpsi] = fem.assema(gr,c,a,frhs);
psi2 = (Kpsi + Mpsi)\Fpsi;


figure(10);
plot(rval,psi1);
title('\psi^{[1]}(r)');
figure(11);
plot(rval,psi2);
title('\psi^{[2]}(r)');


%% Plot B0in animation
gradW = 1;
v1_values = linspace(-1,1,100);

% === POLAR GRID (compute only once) ===
Nr     = 20;
Ntheta = 36;

r_vec     = linspace(0.5,20,Nr);
theta_vec = linspace(0,2*pi,Ntheta+1);
theta_vec(end) = [];

[R,THETA] = meshgrid(r_vec,theta_vec);

X = R.*cos(THETA);
Y = R.*sin(THETA);

% Circle coordinates
theta = linspace(0,2*pi,200);
xc = rcrit*cos(theta);
yc = rcrit*sin(theta);

figure

% ---- First frame ----
v1 = v1_values(1);

fval  = v1*psi1 + gradW*psi2;
frval = diff(fval)./diff(rval);
frval = [0; frval];

rval  = full(rval(:));
fval  = full(fval(:));
frval = full(frval(:));

f    = @(rvar) interp1(rval,fval,rvar,'pchip');
dfdr = @(rvar) interp1(rval,frval,rvar,'pchip');

BR     = -(1./R).*f(R).*cos(THETA);
BTHETA = dfdr(R).*sin(THETA);

BX = BR.*cos(THETA) - BTHETA.*sin(THETA);
BY = BR.*sin(THETA) + BTHETA.*cos(THETA);

q = quiver(X,Y,BX,BY,0.8,'b');
hold on
plot(xc,yc,'r','LineWidth',2)

axis equal
axis([-20 20 -20 20])      % Fix axes for all frames
grid on
xlabel('x')
ylabel('y')

% Uncomment to save a video
% video = VideoWriter('Bfield.mp4','MPEG-4');
% video.FrameRate = 20;
% open(video);

% ---- Animation ----
fig = figure('Position',[100 100 560 420]);
set(fig,'Resize','off');
video = VideoWriter('Bfield.avi','Motion JPEG AVI');
video.FrameRate = 20;
open(video);
for k = 1:length(v1_values)

    v1 = v1_values(k);

    fval  = v1*psi1 + gradW*psi2;
    frval = diff(fval)./diff(rval);
    frval = [0; frval];

    fval  = full(fval(:));
    frval = full(frval(:));

    f    = @(rvar) interp1(rval,fval,rvar,'pchip');
    dfdr = @(rvar) interp1(rval,frval,rvar,'pchip');

    BR     = -(1./R).*f(R).*cos(THETA);
    BTHETA = dfdr(R).*sin(THETA);

    BX = BR.*cos(THETA) - BTHETA.*sin(THETA);
    BY = BR.*sin(THETA) + BTHETA.*cos(THETA);

    % Update quiver data
    q.UData = BX;
    q.VData = BY;

    title(sprintf('Vector field B,  v_1 = %.2f',v1))

    drawnow

    % Save
    writeVideo(video,getframe(gcf));
end

close(video);
%% Parameters
v1 = 1;
gradW = 1;
%% Plot B0in and R0in
% === Plot B-field for first value of v1 ===
%v1 = v1_values(1);

% Compute B0in
fval  = v1*psi1 + gradW*psi2;
frval = diff(fval)./diff(rval);
frval = [0; frval];

fval  = full(fval(:));
frval = full(frval(:));

f    = @(rvar) interp1(rval,fval,rvar,'pchip');
dfdr = @(rvar) interp1(rval,frval,rvar,'pchip');

BR     = -(1./R).*f(R).*cos(THETA);
BTHETA = dfdr(R).*sin(THETA);

BX = BR.*cos(THETA) - BTHETA.*sin(THETA);
BY = BR.*sin(THETA) + BTHETA.*cos(THETA);

% Compute R0in(r)
g1 = H - Hout;
eta=1;
Qin = H.^3/(3*eta); Qout = Hout.^3/(3*eta);
g2 = Qin - Qout;

R0inval  = g1*v1 + g2*gradW;
R0in = @(r) interp1(rval, R0inval, r, 'pchip', 0);   % 0 outside range
RX = R0in(R);
RY = 0.*sin(THETA);

% Plot 
figure;
quiver(X,Y,RX,RY,0.8,'r');
hold on
quiver(X,Y,BX,BY,0.8,'b');
%hold on
% Plot contact line
plot(xc,yc,'r','LineWidth',2)


axis equal; grid on;
axis([-15 15 -15 15])
title(sprintf('{\\color{blue}hatB_0^{in}} and {\\color{red}hatR_0^{in}} for v_1 = %.2f; \\nabla_{\\chi}\\partial_h W_0^{out} = %.2f', v1, gradW))
xlabel('x'); ylabel('y');
%% Plot Q*grad P
QgradPX = -(RX+BX);
QgradPY = -(RY+BY);
figure;
quiver(X,Y,QgradPX,QgradPY,1,'black');
axis equal; grid on;
axis([-15 15 -15 15]);
title(sprintf('- ({\\color{blue}hatB_0^{in}}+{\\color{red}hatR_0^{in}}) for v_1 = %.2f; \\nabla_{\\chi}\\partial_h W_0^{out} = %.2f', v1, gradW))
xlabel('x'); ylabel('y');
%% Plot Stream function Psihat and isolines (level sets)
% === CARTESIAN GRID ===
Nx = 300;
xvec = linspace(-20, 20, Nx);
yvec = linspace(-20, 20, Nx);
[Xc, Yc] = meshgrid(xvec, yvec);

% === POLAR COORDINATES ON GRID ===
Rc     = sqrt(Xc.^2 + Yc.^2);
THETAc = atan2(Yc, Xc);

% === MASK outside r=[0.5, 20] ===
mask = Rc < 0.5 | Rc > 20;

% === EVALUATE cos(theta)*f(r) ===
F_grid = -sin(THETAc) .* f(Rc);   % uses your anonymous function
F_grid(mask) = NaN;                      % hide outside domain

% === PLOT ===
figure;

% Colormap
pcolor(Xc, Yc, F_grid);
shading interp;
%colormap(rdbu);   % or use 'coolwarm', 'rdbu', 'bwr' if no custom cmap
colorbar;
hold on;

% Isolines
contour(Xc, Yc, F_grid, 20, 'k', 'LineWidth', 0.8);

% x axis
plot([-20 20], [0 0], 'black--', 'LineWidth', 1.5);

% Circle at r=rcrit radius of the drop
theta_circ = linspace(0, 2*pi, 500);
plot(rcrit*cos(theta_circ), rcrit*sin(theta_circ), 'r-', 'LineWidth', 1.5);

axis equal; grid on;
xlabel('x'); ylabel('y');
title('-sin(\theta) \cdot f(r)');

% 
% c = x./Qin;
% f=(2-3*min(H)./H);
% Hr = diff(H)./diff(x);
% Hr=[Hr; Hr(end)];
% frhs = -(x.*f.*Hr)./Qin;
% fem=p.pdeo.fem;
% gr=p.pdeo.grid; 
% [K_psi,M_psi,F_psi] = fem.assema(gr,c,0,frhs);
% psi1 = K_psi \ F_psi;



% Compute residual without FEM
% psix  = gradient(psi, x);        % first derivative
% psixx = gradient(psix, x);       % second derivative
% res = (x.^2).*psixx + x.*psix.*(1-x.*drlogQin) - psi + (x.^2).*gradient(h,x);
% norm(res)

%p=cont(p,10); 
%branch = p.branch;
%writematrix(branch,'1D_c0_continuation_22_04_26.txt');