%% demo for pBC 1D, clear workspace 
close all; keep pphome; 
%% cell 1: init
p=[]; par=[0.5 0 0]; % initial concentration, chemical potential (lagrange mult mass)
p=chinit(p,15,2500,par); p.sw.qjac=0; %p.sw.verb=2;
%% Continuation parameters
p.nc.ilam = [2, 3];
p.nc.nq=1;
p.nc.lammax=50; p.sol.ds=0.1; p.nc.dsmax=0.3;
%% First branch continuation
p=setfn(p,'tr'); p=findbif(p,100);
%% Switch to droplet branch
p.nc.lammax=100;
p.sol.ds=+0.1;
p.nc.dsmax=0.1;
p=swibra('tr','bpt1','b1',+0.01); p=cont(p,20); 
p.sw.bifcheck=0;
p.sw.foldcheck=1;
%p.sol.ds=-0.01;
%p=cont(p,100);
%% Save observables

branch = full(p.branch);
writematrix(branch,'2D_continuation.txt');
Hmax=branch(7,:);
Hout=branch(8,:);
Hout_ = branch(13,:);
Psi1lim = branch(14,:);
c0 = branch(15,:);
r0 = branch(16,:);
I1=branch(9,:);I2=branch(10,:);
figure(20)
plot(Hout_,I1);hold on;plot(Hout_,I2); xlabel('H_0^{out}');legend('I_1','I_2');ylim([0,4]);grid();
I1=branch(11,:);I2=branch(12,:);
figure(21)
plot(Hout_,I1);hold on;plot(Hout_,I2); xlabel('H_0^{out}');legend('I_1','I_2');ylim([0,4]);grid();
% figure(80)
% hold on
% plot(Hout_, Hmax);
% xlabel('H_0^{out}');
% ylabel('H_{max}');
% grid("on");
% xlim([1.05, 1.3]);
figure(81);
plot(Hout_, Hout);
hold on
plot(Hout_, Hout_,'--');
xlabel('H_0^{out} measured');
ylabel('H_0^{out} formula');
grid("on");
xlim([1.01, 1.3]);
ylim([1.01, 1.3]);

figure(91);
plot(Hout_, Psi1lim);
hold on
ylabel('\psi^{[1]}(\infty)');
xlabel('H_0^{out} measured');
figure(92);
plot(Hout_, r0);
hold on
ylabel('drop radius');
xlabel('H_0^{out} measured');
%%
H = p.u(1:p.nu); % Exclude the parameters of the pde
Hout = Hout_;
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
%eta=1;
%Qin = H.^3/(3*eta);

c = (H.^3).^(-1);
a = 1 ./ ((1e-12+rval.^2) .* (H.^3));    % 1e-12 To eliminate the divergence at r=0

%fac  = 1 - 3*g./H;
Hr = diff(H)./diff(rval);
Hr = [0; Hr];
%frhs = -(fac .* Hr) ./ Qin;
G = (3*Hout./H-2);
frhs = -G.*(Hr./(H.^3));

fem=p.pdeo.fem;
gr=p.pdeo.grid;
[Kpsi,Mpsi,Fpsi] = fem.assema(gr,c,a,frhs);
psi1 = (Kpsi + Mpsi)\Fpsi;
% Dirichlet right at x=L
A = Kpsi + Mpsi;
A(end,:) = 0;
A(:,end) = 0;
A(end, end) = 1;
Fpsi(end) = 0;
%psi1 = A\Fpsi;

% Psi 2

%Qin  = H.^3/(3*eta);
%Qout = Hout.^3/(3*eta);
% c(r) and a(r) are the same! Only frhs is different

%frhs = -3 .* (Hout.^3 ./ H.^4) .* Hr;
G = (Hout.^3)./H;
frhs = -G.*(Hr./(H.^3));

[Kpsi,Mpsi,Fpsi] = fem.assema(gr,c,a,frhs);
psi2 = (Kpsi + Mpsi)\Fpsi;
% Dirichlet right at x=L
A = Kpsi + Mpsi;
A(end,:) = 0;
A(:,end) = 0;
A(end, end) = 1;
Fpsi(end) = 0;
%psi2 = A\Fpsi;


figure(10);
plot(rval,psi1);
title('\psi^{[1]}(r)');
figure(11);
plot(rval,psi2);
title('\psi^{[2]}(r)');

psi1(end)
psi2(end)
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
gradW = 1;  % This is gradw*ha2*eta^-1
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
R0inval  = (H-Hout)*v1 + (gradW/3)*(H.^3-Hout.^3);
R0in = @(r) interp1(rval, R0inval, r, 'pchip', 0);   % 0 outside range
RX = R0in(R);
RY = 0.*sin(THETA);

% Plot 
% Same scale for the two vector fields
%mag1 = hypot(RX,RY);
%mag2 = hypot(BX,BX);
%maxMag = max([mag1(:); mag2(:)]);

figure;
R0 = hypot(RX,RY);
pcolor(X,Y,log10(R0))
shading interp
alpha(0.3)
hold on

quiver(X,Y,RX,RY,1,'r')

axis equal tight
colormap(turbo)
cb = colorbar;
cb.Label.String = 'log_{10}(|hatR_0^{in}|)';
hold on;
plot(xc,yc,'black','LineWidth',2)
axis equal; grid on;
axis([-15 15 -15 15])
title(sprintf('{{\\color{red}h_a^{-1}hatR_0^{in}} }for v_1 = %.2f; (h_a^2\\eta^{-1}\\nabla_{\\chi}\\partial_h W_0^{out}) = %.2f', v1, gradW))
xlabel('x'); ylabel('y');



figure;
B0 = hypot(BX,BY);
pcolor(X,Y,log10(B0))
shading interp
alpha(0.3)
hold on
quiver(X,Y,BX,BY,0,'b');
colormap(turbo)
cb = colorbar;
cb.Label.String = 'log_{10}(|hatB_0^{in}|)';
hold on
% Plot contact line
plot(xc,yc,'black','LineWidth',2)


axis equal; grid on;
axis([-15 15 -15 15])
title(sprintf('{\\color{blue}h_a^{-1}hatB_0^{in}} for v_1 = %.2f; (h_a^2\\eta^{-1}\\nabla_{\\chi}\\partial_h W_0^{out}) = %.2f', v1, gradW))
xlabel('x'); ylabel('y');
%% Plot Q*grad P
QgradPX = (BX+RX);
QgradPY = (BY+RY);
figure;
quiver(X,Y,QgradPX,QgradPY,1,'black');
hold on
% Plot contact line
plot(xc,yc,'black','LineWidth',2)

axis equal; grid on;
axis([-15 15 -15 15]);
title(sprintf(' ({\\color{blue}h_a^{-1}hatB_0^{in}}+{\\color{red}h_a^{-1}hatR_0^{in}}) for v_1 = %.2f; (h_a^2\\eta^{-1}\\nabla_{\\chi}\\partial_h W_0^{out}) = %.2f', v1, gradW))
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