function [nep] = nep_wg_generator(nx, nz, options)
%NEP_WG_GENERATOR Generate the waveguide eigenvalue problem
%[nep] = nep_wg_generator(nx, nz, options)
%
%  This is an optimized version of NEP_WG_GENERATOR which generates the
%  action of the matrices implicitly without constructing them.
%
%  Input:
%    nx = number of points in x-direction
%    nz = number of points in z-direction
%    options = struct of different options to choose from. Valid options:
%      delta = extra space added to x_minus and x_plus for the cut (default is 0).
%      wg = tells you what wave gudie it is.  TAUSCH (default)
%                                             CHALLENGE
%
%  Output:
%    nep = a struct holding a variety of functions for generating matrices
%    related to the NEP
%
%  EXMAPLE EXECUTION
% % n = 100;
% % nep_options.delta = 0.1;
% % nep_options.wg = 'TAUSCH';
% % nep = nep_wg_generator(n, n, nep_options);

% INPUT AND OPTION HANDLING
options_are_given = exist('options', 'var');
if(options_are_given && isfield(options, 'delta'))
  delta = options.delta;
else
  delta = 0;
end
if(options_are_given && isfield(options, 'wg'))
  wg = options.wg;
else
  wg = 'TAUSCH';
end



nep.nx = nx;
nep.nz = nz;

% SET DOMAIN SIZE
if(strcmp(wg, 'TAUSCH'))
  xm = 0;
  xp = 2/pi + 0.4;
elseif(strcmp(wg, 'CHALLENGE'))
  xm = -1;
  xp = 1;
else
  error('Not a vaild Wave Guide');
end
xm = xm - delta;   nep.xm = xm;
xp = xp + delta;   nep.xp = xp;

% GENERATION OF THE DOMAIN DISCRETIZATION
% domain (First generate including the boundary)
X = linspace(xm, xp, nx+2);
Z = linspace(0, 1, nz+1);
% Removing the boundary
X = X(2:end-1);
Z = Z(2:end);

% DISCRETIZATION STEP
hx = X(2)-X(1);
hz = Z(2)-Z(1);
nep.hx = hx;
nep.hz = hz;

ex = ones(nx,1);
ez = ones(nz,1);

% DISCRETIZATION OF THE SECOND DERIVATIVE
Dxx = spdiags([ex -2*ex ex], -1:1, nx, nx); 
Dzz = spdiags([ez -2*ez ez], -1:1, nz, nz); 

% IMPOSE PERIODICITY
Dzz(1, end) = 1;
Dzz(end, 1) = 1;

Dxx = Dxx/(hx^2);
Dzz = Dzz/(hz^2);

% DISCRETIZATION OF THE FIRST DERIVATIVE
Dz  = spdiags([-ez ez], [-1 1], nz, nz);

% IMPOSE PERIODICITY
Dz(1,end) = -1;
Dz(end,1) = 1;

Dz = Dz/(2*hz);

% WAVENUMBER
if(strcmp(wg, 'TAUSCH'))
  k1 = sqrt(2.3)*pi;
  k2 = sqrt(3)*pi;
  k3 = pi;

  k = @(x,z) ...
          k1*(x <= 0) +                           ...
          k2*(x>0).*(x<=2/pi) +                   ...
          k2*(x>2/pi).*(x<=2/pi+0.4).*(z>0.5)+    ...
          k3*(x>2/pi).*(z<=0.5).*(x<=2/pi+0.4)+   ...
          k3*(x>2/pi+0.4);
elseif(strcmp(wg, 'CHALLENGE'))
  k1 = sqrt(2.3)*pi;
  k2 = 2*sqrt(3)*pi;
  k3 = 4*sqrt(3)*pi;
  k4 = pi;
  
  k = @(x,z) ...
          k1 *(x<=-1) +                                              ...
          k4 *(x>1) +                                                 ...
          k4 *(x>(-1+1.5)) .* (x<=1) .* (z<=0.4) +                     ...
          k3 *(x>(-1+1)) .* (x<=(-1+1.5)) +                             ...
          k3 *(x>(-1+1.5)) .* (x<=1) .* (z>0.4) +                      ...
          k3 *(x>-1) .* (x<=(-1+1)) .* (z>0.5) .* (z-x/2<=1) +        ...
          k2 *(x>-1) .* (x<=(-1+1)) .* (z>0.5) .* (z-x/2>1) +         ...
          k3 *(x>-1) .* (x<=(-1+1)) .* (z<=0.5) .* (z+x/2>0) +        ...
          k2 *(x>-1) .* (x<=(-1+1)) .* (z<=0.5) .* (z+x/2<=0);
else
  error('Not a vaild Wave Guide');
end

[x_mesh, z_mesh] = meshgrid(X, Z);
nep.K = k(x_mesh, z_mesh).^2;


% BUILD THE FIRST BLOCK Q
Ix = speye(nx);   Iz=speye(nz);
nep.Q = @(gamma,x) reshape((Dzz + 2*gamma*Dz + gamma^2*Iz)*reshape(x,nz,nx) + reshape(x,nz,nx)*Dxx + nep.K.*reshape(x,nz,nx), nx*nz, 1);

% BUILD THE DERIVATIVE OF Q
nep.p_Q = @(gamma,x) reshape((2*Dz + 2*gamma*Iz)*reshape(x,nz,nx), nx*nz, 1);



% BUILD THE SECOND BLOCK C1
e1 = Ix(:,1);
en = Ix(:,end); 
nep.C1 = [kron(e1,Iz) kron(en,Iz)]/(hx^2);


% BUILD THE THIRD BLOCK C2^T
d1 = 2/hx;              nep.d1=d1;
d2 = -1/(2*hx);         nep.d2=d2;
vm = sparse(1,nx);
vm(1) = d1;
vm(2) = d2;
vp = sparse(1,nx);
vp(end) = d1;
vp(end-1) = d2;

nep.C2T = [kron(vm,Iz); kron(vp,Iz)];


% FOURIES OPERATORS FOR MATRICES
p = (nz-1)/2;
bb = exp(-2i*pi*((1:nz)-1)*(-p)/nz).';  % scaling to do after FFT
nep.R = @(X) flipud((bb*ones(size(X,2),1)').*fft(X));
bbinv = 1./bb; % scaling to do before inverse FFT
nep.Rinv = @(X) ifft((bbinv*ones(size(X,2),1)').*flipud(X));


% FUNCTIONS (AND DERIVATIVES) RELATED TO THE DtN-MAP
Km = k(-Inf, 1/2);
Kp = k(Inf, 1/2);
d0 = -3/(2*hx);         nep.d0=d0;
a = ones(nz,1);
b = 4*pi*1i*(-p:p).';
cM = Km^2-4*pi^2*((-p:p).^2).';
cP = Kp^2-4*pi^2*((-p:p).^2).';

betaM = @(gamma) a*gamma^2+b*gamma+cM;
betaP = @(gamma) a*gamma^2+b*gamma+cP;

signM = 1i*sign(imag(betaM(-1-1i)));                % OBS! LEFT HALF-PLANE!
signP = 1i*sign(imag(betaP(-1-1i)));                % OBS! LEFT HALF-PLANE!

nep.sM = @(gamma) signM.*sqrt(betaM(gamma))+d0;
nep.sP = @(gamma) signP.*sqrt(betaP(gamma))+d0;

nep.p_sM = @(gamma) signM.*(2*a*gamma+b)./(2*sqrt(a*gamma^2+b*gamma+cM));
nep.p_sP = @(gamma) signP.*(2*a*gamma+b)./(2*sqrt(a*gamma^2+b*gamma+cP));

% BUILD THE FOURTH BLOCK P
nep.P=@(gamma,x)  [nep.R(nep.Rinv(x(1:end/2)).*nep.sM(gamma))      ;
                   nep.R(nep.Rinv(x(end/2+1:end)).*nep.sP(gamma))  ];

% BUILD THE DERIVATIVE OF P
nep.p_P=@(gamma,x)  [nep.R(nep.Rinv(x(1:end/2)).*nep.p_sM(gamma))      ;
                     nep.R(nep.Rinv(x(end/2+1:end)).*nep.p_sP(gamma))  ];


% ASSEMBLE THE MATRIX M
nep.M = @(gamma,x) [nep.Q(gamma,x(1:nx*nz)) + nep.C1*x(nx*nz+1:nx*nz+2*nz)        ;          
                    nep.C2T*x(1:nx*nz)      + nep.P(gamma, x(nx*nz+1:nx*nz+2*nz) )  ];
          
% ASSEMBLE DERIVATIVE OF THE MATRIX M
nep.p_M = @(gamma,x) [nep.p_Q(gamma,x(1:nx*nz));          
                      nep.p_P(gamma, x(nx*nz+1:nx*nz+2*nz))    ];


% BUILD THE INVERSE OF THE MATRIX P
nep.Pinv=@(gamma,x)  [  nep.R(nep.Rinv(x(1:end/2))./nep.sM(gamma))      ;
                        nep.R(nep.Rinv(x(end/2+1:end))./nep.sP(gamma))  ];

nep.Pm_inv=@(gamma, x) nep.R(nep.Rinv(x)./nep.sM(gamma));                     
nep.Pp_inv=@(gamma, x) nep.R(nep.Rinv(x)./nep.sP(gamma));                     


%ASSEMBLE THE SCHUR COMPLEMENT
nep.S = @(gamma,x) nep.Q(gamma,x) - nep.C1*nep.Pinv(gamma,nep.C2T*x);

fprintf('nep_wg_generator is done!\n')
end