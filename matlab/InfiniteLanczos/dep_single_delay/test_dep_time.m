close all
clear all
clc

rng(1)


%fftw('dwisdom',[]);
fftw('planner','measure');
%fftw('planner','estimate');


n=1e4;  % size of the problem
m=100;  % number of iterations of Infinite Lanczos

%m=300;


% % generate the problem
A0 = n^2*spdiags(rand(n,3), -1:1, n, n);    A0=A0+A0';
A1 = n^2*spdiags(rand(n,3), -1:1, n, n);    A1=A1+A1';


nep = load_dep( A0, A1 ); nep.Md =@(j) nep.Md(j);
nep.Md_lin_comb=@(X,j) -X(:,1)+A1*(sum(bsxfun(@times, X(:,1:j),(-1).^(1:j)),2));

[alpha, beta] = G(m+1);
nep.alpha=alpha;
nep.beta=beta;


v1=rand(n,1);   

% coefficients matrix (factorials)
cc=gen_coeffs(2*(m+1)); cc=cc(1:m+1,1:m+1); nep.cc=cc;

t1=cputime; [ T, Omega, W, H ] = Infinite_Lanczos_reduced( v1, m, nep ); t1=cputime-t1;

% Generate the error history by extracting the eig approx from the 
% projection matrix (as expected does not work well).
%[ ~, D, err ] = error_hyst( V, T, H, Omega, nep );

% Solve the same problem with TIAR
nep.Md=@(j,v) nep.Md(j)*v; % adapt the interface for TIAR
opts.maxit=m; opts.tol=1e-12; opts.disp=1;  k=Inf;
opts.sigma=0;    opts.gamma=1;  opts.p=Inf;
t2=cputime; [ ~, E, ~ ] = tiar( nep, k, opts ); t2=cputime-t2;

fprintf("INFLAN %f seconds\n",t1)
fprintf("TIAR %f seconds\n",t2)