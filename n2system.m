% continuum parameters
lamy = @(x,y) ones(size(x));
muy1 = @(x) 2*ones(size(x));
muy2 = @(x) 1*ones(size(x));
sigy = @(x,y,h) x.*(x+1).*(y-0.5).*x.^2.*(h-0.5);
Wy1 = @(x,y) x.*(x+1).*(y-0.5).*exp(x);
Wy2 = @(x,y) x.*(x+1).*(y-0.5).*exp(x);
thy1 = @(x,y) -3*y.*(y-1).*ones(size(x));
thy2 = @(x,y) -2*y.*(y-1).*ones(size(x));
qy1 = @(y) 8*(y-0.5);
qy2 = @(y) -8*(y-2);
psi12y = @(x) 0*ones(size(x));
psi21y = @(x) 0*ones(size(x));
% continuum kernels solutions and control law
L22c = @(x,xi) -2*exp(2*(x-xi));
L12p2c = @(x,xi) -2*exp(x-2*xi); % p=2, L12p1 = 0;
K1p1c = @(x,xi,y) -thy1(x,y)/3; % p=1
K1p2c = @(x,xi,y) -exp(x-2*xi).*thy1(x,y)/3; % p=2
K2c = @(x,xi,y) -exp(2*(x-xi)).*thy2(x,y)/2;

% discretization
xn = 257; % grid points for x
xg = 256; % grid points for u and v (excluding boundary conditions)
x = linspace(0,1,xn); % grid for x
x0 = x(1:xg); x1 = x(2:xg+1); % grids for u and
x0h = x0(1:xg/2); % first half of x0
x1h = x1(1:xg/2); xh1 = x1(xg/2+1:xg); % first and second halfs of x1

% D is backward difference, -D' is forward difference
D = eye(xg) - diag(ones(1,xg-1), -1);

% considered values of n
nn = [2 6 10 10];
% cell arrays for storing the controls
U1c = cell(1,4); U2c = U1c;

for ni = 1:4
n = nn(ni); y = linspace(1/n,1,n); % n and grid for y
% n+m parameters parameters (sample from continuum for given y)
lam = @(i,x) lamy(x,y(i));
mu1 = @(x) muy1(x);
mu2 = @(x) muy2(x);
sig = @(i,j,x) sigy(x,y(i),y(j));
W1 = @(i,x) Wy1(x,y(i));
W2 = @(i,x) Wy2(x,y(i));
th1 = @(j,x) thy1(x,y(j));
th2 = @(j,x) thy2(x,y(j));
psi12 = psi12y(x0);
psi21 = psi21y(x0);

q1 = qy1(y); q2 = qy2(y);

% finite difference approximation into \dot x = Ax + BU, where
A = zeros((n+2)*xg); % system operator (initialization)
B = zeros((n+2)*xg, 2);
B((n+1)*xg,1) = mu1(1)*xg; % control operator v1(1) = U1
B((n+2)*xg,2) = mu2(1)*xg; % control operator v2(1) = U2

% fill A
for k = 1:n
  Ik = (k-1)*xg+1:k*xg; % index set
  A(Ik,Ik) = -xg*lam(k,x1).*D; % tranport term
  A(Ik(1),n*xg+1) = xg*lam(k,0)*q1(k); % boundary condition
  A(Ik(1),(n+1)*xg+1) = xg*lam(k,0)*q2(k); % boundary condition
  for l = 1:n
    Il = (l-1)*xg+1:l*xg; % index set
    A(Ik,Il) = A(Ik,Il) + diag(sig(k,l,x1))/n; % sigma terms (scale 1/n)
  end
  Il = n*xg+1:(n+1)*xg; % index set
  A(Ik,Il) = A(Ik,Il) + diag(W1(k,x0)); % W1 term
  Il = (n+1)*xg+1:(n+2)*xg; % index set
  A(Ik,Il) = A(Ik,Il) + diag(W2(k,x0)); % W2 term
end
k=n+1;
Ik = (k-1)*xg+1:k*xg; % index set
A(Ik,Ik) = -xg*mu1(x0).*D'; % transport term
for l = 1:n
  Il = (l-1)*xg+1:l*xg; % index set
  A(Ik,Il) = A(Ik,Il) + diag(th1(l,x1))/n; % theta1 terms (scale 1/n)
end
k=n+2;
Ik = (k-1)*xg+1:k*xg; % index set
A(Ik,Ik) = -xg*mu2(x0).*D'; % transport term
for l = 1:n
  Il = (l-1)*xg+1:l*xg; % index set
  A(Ik,Il) = A(Ik,Il) + diag(th2(l,x1))/n; % theta2 terms (scale 1/n)
end
Il = xg*n+1:xg*(n+1);
A(Ik,Il) = diag(psi21);
A(Il,Ik) = diag(psi12);
% max(real(eig(A))) % open-loop stability check

% control law
if ni == 3 % use 10+m kernels solved in n2kernelsolver.m
  run n2kernelsolver.m
  U1 = @(z) (trapz(K1n.*z(1:n*xg))/n + ...
    + trapz([L11n; L12n].*z(n*xg+1:(n+2)*xg)))/xg;
  U2 = @(z) (trapz(K2n.*z(1:n*xg))/n + ...
    trapz([L21n; L22n].*z(n*xg+1:(n+2)*xg)))/xg;
else % % compute control gains based on continuum kernels
  K1p1n = zeros(xg*n, 1); K1p2n = K1p1n; K2n = K1p1n;
  for k = 1:n
    K1p1n((k-1)*xg+xg/2+1:k*xg) = K1p1c(1, xh1, y(k))';
    K1p2n((k-1)*xg+1:(k-1)*xg+xg/2) = K1p2c(1, x1h, y(k))';
    K2n((k-1)*xg+1:k*xg) = K2c(1, x1, y(k))';
  end
  U1 = @(z) (trapz(K1p1n.*z(1:n*xg))/n + ...
    trapz(K1p2n.*z(1:n*xg))/n + ...
    trapz(L12p2c(1,x0h').*z((n+1)*xg+1:(n+1)*xg+xg/2)))/xg;
  U2 = @(z) (trapz(K2n.*z(1:n*xg))/n + ...
    trapz(L22c(1,x0').*z((n+1)*xg+1:(n+2)*xg)))/xg;
end

% simulate
opts = odeset('Jacobian', A); % should expedite computations
T = 5; % simulation end time
z0 = [kron(q1'+q2', ones(xg,1)); ones(2*xg,1)]; % initial condition
% solve; ode45 (ok fast), ode23 (faster, accuracy?), others much slower
sol = ode45(@(t,z) A*z + B*[U1(z); U2(z)], [0, T], z0, opts);
tg = 513; 
TT = linspace(0,T,tg); % time grid (for plotting)
usol = deval(sol, TT); % evaluate solution at tg
% compute and store controls
Usol1 = zeros(1,tg); % evaluate input at tg
Usol2 = Usol1;
for k=1:tg
  Usol1(k) = U1(usol(:,k));
  Usol2(k) = U2(usol(:,k));
end
U1c{ni} = Usol1; U2c{ni} = Usol2;
end

%% plot results
figure(1) % inputs
subplot(211)
plot(TT, U1c{1},'linewidth',2,'DisplayName', '$U_2^1(t)$')
hold on
plot(TT, U1c{2},'linewidth',2,'DisplayName', '$U_6^1(t)$')
plot(TT, U1c{3},'linewidth',2,'DisplayName', '$U_{10}^1(t)$')
plot(TT, U1c{4},'linewidth',2,'DisplayName', '$U^1_{\mathrm e}(t)$')
hold off
set(gca,'ytick',-15:5:5,'xtick',1:6,'xticklabel','{}')
set(gca,'tickdir', 'out', 'fontsize',12)
set(gca,'fontsize',11)
ylabel('$U^1(t)$', 'interpreter','latex', 'fontsize', 12, ...
  'rotation', 0)
legend('interpreter', 'latex','fontsize',12,'location','southeast', ...
  'numcolumns',2)
ylim([-12 2])
xlim([0 5])
set(gca,'position',get(gca,'position')+[0 0.02 0.06 0.02])
subplot(212)
plot(TT, U2c{1},'linewidth',2, 'DisplayName', '$U_2^2(t)$')
hold on
plot(TT, U2c{2},'linewidth',2, 'DisplayName', '$U_6^2(t)$')
plot(TT, U2c{3},'linewidth',2, 'DisplayName', '$U_{10}^2(t)$')
plot(TT, U2c{4},'linewidth',2, 'DisplayName', '$U^2_{\mathrm e}(t)$')
hold off
set(gca,'ytick',-15:5:5)
set(gca,'tickdir', 'out', 'fontsize',12)
set(gca,'fontsize',11)
xlabel('$t$', 'interpreter', 'latex', 'fontsize',12)
ylabel('$U^2(t)$', 'interpreter','latex', 'fontsize', 12, ...
  'rotation', 0)
set(gca,'position',get(gca,'position')+[0 0.06 0.06 0.02])
legend('interpreter', 'latex','fontsize',12,'location','southeast', ...
  'numcolumns',2)
ylim([-13 4])
xlim([0 5])

%%
figure(2) % states
pn1 = n;
pn2 = n+1;
subplot(1,2,1)
surf(TT, x1, usol((pn1-1)*xg+1:pn1*xg,:))
shading interp
xlabel('$t$','interpreter', 'latex','fontsize', 12)
ylabel('$x$','interpreter', 'latex','fontsize', 12)
set(gca,'ytick',[0 1],'fontsize',11)
set(gca,'xtick',[0 5],'fontsize',11)
xlim([0 5])
zlabel(['$u^{',num2str(pn1),'}(t,x)$'],'interpreter','latex','fontsize',12)
subplot(1,2,2)
surf(TT, x1, usol((pn2-1)*xg+1:pn2*xg,:))
shading interp
xlabel('$t$','interpreter', 'latex','fontsize', 12)
ylabel('$x$','interpreter', 'latex','fontsize', 12)
set(gca,'ytick',[0 1], 'fontsize',11)
set(gca,'xtick',[0 5], 'fontsize',11)
xlim([0 5])
zlabel('$v^{1}(t,x)$','interpreter', 'latex','fontsize', 12)
