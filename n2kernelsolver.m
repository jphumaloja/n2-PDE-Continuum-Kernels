%% precompute everything, initialize
n = 10; y = linspace(1/n,1,n); % grid for y
xn = 257; x = linspace(0,1,xn);
% continuum paramters
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
% aux, for precomputing, reduces function calls
q1 = qy1(y);
q2 = qy2(y);
S = zeros(n,n);
WM1 = zeros(n,xn);
THM1 = WM1; WM2 = WM1; THM2 = WM1;
L3 = zeros(n,n,xn); U3 = L3; L2 = L3; U2 = U3;
for k=1:xn
  for s=1:n
    S(s,:) = sigy(x(k),y(1:n),y(s)); % S
  end
  % precompute thy and Wy
  THM1(:,k) = thy1(x(k),y(1:n));
  WM1(:,k) = Wy1(x(k),y(1:n));
  THM2(:,k) = thy2(x(k),y(1:n));
  WM2(:,k) = Wy2(x(k),y(1:n));
  % precompute LU-decompositions of (2+S) and (3+S)
  [L2(:,:,k), U2(:,:,k)] = lu(2*eye(n)-S/(n*xn));
  [L3(:,:,k), U3(:,:,k)] = lu(3*eye(n)-S/(n*xn));
end
kd1 = zeros(1,1,n); kd2 = zeros(1,1,n);
% kernels K and L
K1 = zeros(xn,xn,n); % 3D matrix, components x, xi, i
K2 = K1;
% L as 2D matrices, compoents x, xi (initial guess to zero)
L22 = zeros(xn,xn); L21 = L22; L11 = L22;
L12p1 = L22; L12p2 = L22; % only L12 is split (hopefully this works)
% boundary condition for K on xi = x
for k = 1:xn % probably some better way to do this (3D matrices suck)
  for kn = 1:n
    K1(k,k,kn) = -thy1(x(k),y(kn))/3;
    K2(k,k,kn) = -thy2(x(k),y(kn))/2;
  end
end
%%
% solve iteratively until changes are smaller than the chosen tolerance
iter = 400; % max iter
tol = 1e-7; % cutout 
dn = zeros(1,iter);
% k=x, m=xi
for mm = 1:iter
  % disp(mm) % to show the iteration 
  ok1 = K1; ok2 = K2;
  ol11 = L11; ol21 = L21; ol22 = L22;
  ol12p1 = L12p1; ol12p2 = L12p2;
  % solve K based on L
  for k = 2:xn
    for m = 1:k-1
      % compute everything in 2D, then tranform to 3D
      % crude finite-difference approximation of the kernel equations
      kv1 = 2*K1(k-1,m,:) + K1(k,m+1,:);
      kv1 = kv1(:);
      if m < k/2 % characteristic line, might need some more work
        kd1(1,1,:) = U3(:,:,k)\(L3(:,:,k)\ ...
          (kv1 + (THM1(:,m)*L11(k,m) + THM2(:,m)*L12p2(k,m))/xn));
      else
        kd1(1,1,:) = U3(:,:,k)\(L3(:,:,k)\ ...
          (kv1 + (THM1(:,m)*L11(k,m) + THM2(:,m)*L12p1(k,m))/xn));
      end
      K1(k,m,:) = kd1;
      kv2 = K2(k-1,m,:) + K2(k,m+1,:);
      kv2 = kv2(:);
      kd2(1,1,:) = U2(:,:,k)\(L2(:,:,k)\ ...
          (kv2 + (THM1(:,m)*L21(k,m) + THM2(:,m)*L22(k,m))/xn));
      K2(k,m,:) = kd2;
    end
  end  
  % solve L based on K
  L11(:,1) = 0.5*K1(:,1,1)*q1(1);
  L12p2(:,1) = K1(:,1,1)*q2(1);
  L22(:,1) = K2(:,1,1)*q2(1);
  for k = 2:n
    L11(:,1) = L11(:,1) + 0.5*K1(:,1,k)*q1(k);
    L12p2(:,1) = L12p2(:,1) + K1(:,1,k)*q2(k);
    L22(:,1) = L22(:,1) + K2(:,1,k)*q2(k);
  end
  L11(:,1) = L11(:,1)/n;
  L12p2(:,1) = L12p2(:,1)/n;
  L22(:,1) = L22(:,1)/n;
  for k = 2:xn
    for m = 2:k-1
      % crude finite-difference approximation of the kernel equation
      k1v = K1(k,m,:); k2v = K2(k,m,:);
      % may need to account for the characteristic line with L11
      L11(k,m) = 0.25*(2*L11(k-1,m) + 2*L11(k,m-1) ...
        + sum(WM1(:,m).*kv1(:))/(n*xn));
      L22(k,m) = 0.5*(L22(k-1,m) + L22(k,m-1) ...
        + sum(WM2(:,m).*kv2(:))/(n*xn));
      if m <= k/2 % L12p2 only computed up to xi = 0.5x
        L12p2(k,m) = 1/3*(2*L12p2(k-1,m) + L12p2(k,m-1) ...
        + sum(WM1(:,m).*kv2(:))/(n*xn));
      end
    end
    % at the boundary xi=x, use m=k-1 as central point
    kv1 = K1(k,k-1,:); kv2 = K2(k,k-1,:);
    L11(k,k) = 2*L11(k,k-1)-L11(k-1,k-1)+0.5*sum(WM1(:,k-1).*kv1(:)/(n*xn));
    L22(k,k) = 2*L22(k,k-1)-L22(k-1,k-1)+sum(WM2(:,k-1).*kv2(:)/(n*xn));
  end
  % remaining L21 and L12p1 (starting from xi = x with 0 boundary condtion)
  for k = 2:xn
    for m = 1:k-1
      % crude finite-difference approximation of the kernel equations
      k1v = K1(k,m,:); k2v = K2(k,m,:);
      L21(k,m) = 1/3*(L21(k-1,m) - 2*L21(k,m+1) ...
        + sum(WM2(:,m).*kv1(:))/(n*xn));
      if m >= k/2 % L12p1 only computed down to xi = 0.5x
        L12p1(k,m) = 1/3*(2*L12p1(k-1,m) - L12p1(k,m+1) ...
          + sum(WM1(:,m).*kv2(:))/(n*xn));
      end
    end
  end  
  % check changes, terminate if small, otherwise repeat
  dn(mm) = norm(K1(:)-ok1(:),inf) + norm(K2(:)-ok2(:),inf) ...
    + norm(L11(:)-ol11(:),inf) + norm(L12p1(:)-ol12p1(:),inf) ...
    + norm(L12p2(:)-ol12p2(:),inf) + norm(L21(:)-ol21(:),inf) ...
    + norm(L22(:)-ol22(:),inf);
  if dn(mm) < tol
    % show terminating iteration
    disp(['kernelsolver converged at iteration ',num2str(mm)]) 
    break
  end
end
% evaluate control gains (at x = 1)
L11n = L11(xn,1:xn-1)';
L12n = [L12p2(xn,1:(xn-1)/2)'; L12p1(xn,(xn-1)/2+1:xn-1)'];
L21n = L21(xn,1:xn-1)'; L22n = L22(xn,1:xn-1)';
K1n = zeros(n*(xn-1),1); K2n = K1n;
for k = 1:n
  K1n((k-1)*(xn-1)+1:k*(xn-1)) = K1(xn,2:xn,k)';
  K2n((k-1)*(xn-1)+1:k*(xn-1)) = K2(xn,2:xn,k)';
end