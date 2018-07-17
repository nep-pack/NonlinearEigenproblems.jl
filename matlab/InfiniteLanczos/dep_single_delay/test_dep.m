close all
clear all
clc

rng(1)

n=200;  % size of the problem
m=100;  % number of iterations of Infinite Lanczos

% generate the problem
A0 = n^2*spdiags(rand(n,3), -1:1, n, n);    A0=A0+A0';
A1 = n^2*spdiags(rand(n,3), -1:1, n, n);    A1=A1+A1';

nep.Md_lin_comb=@(X,j) -X(:,1)+A1*(sum(bsxfun(@times, X(:,1:j),(-1).^(1:j)),2));
nep.n=n; nep.A1=A1; 
[L,U,P,Q] = lu(A0+A1);
nep.M_inv=@(v) Q*(U\(L\(P*v)));
nep.M =@(l,v) -l*v+A0*v+exp(-l)*(A1*v);
nep.resnorm=@(l,v) norm(nep.M(l,v))/((abs(l)+norm(A0,inf)+abs(exp(-l))*norm(A1,inf))*norm(v));



v1=rand(n,1);   

% coefficients matrix (factorials)
cc=gen_coeffs(2*(m+1)); cc=cc(1:m+1,1:m+1); nep.cc=cc;
t1=cputime; [ T, Omega, W, H ] = Infinite_Lanczos_reduced( v1, m, nep ); t1=cputime-t1;

% Generate the error history by extracting the eig approx from the 
% projection matrix (as expected does not work well).
%[ ~, D, err ] = error_hyst( V, T, H, Omega, nep );

% Solve the same problem with TIAR
opts.maxit=m; opts.tol=1e-12; opts.disp=1;  k=Inf;
opts.sigma=0;    opts.gamma=1;  opts.p=1;
t2=cputime; [ ~, E, ~ ] = tiar( nep, k, opts ); t2=cputime-t2;

% plot the eigenvalues computed with TIAR
figure(2)
plot(real(E),imag(E),'*'); hold on
theta=linspace(0,2*pi,100);
plot(cos(theta),sin(theta),'k-'); % plot the unit disk


% extract the first block row from Infinite Lanczos.  
VV=W; %[VV,~]=qr(VV,0);
% construct the projected problem
A0=VV'*(A0*VV); A1=VV'*(A1*VV);
nep.Md_lin_comb=@(X,j) -X(:,1)+A1*(sum(bsxfun(@times, X(:,1:j),(-1).^(1:j)),2));
nep.n=length(A0); nep.A1=A1; 
[L,U] = lu(A0+A1);
nep.M_inv=@(v) (U\(L\v));
nep.M =@(l,v) -l*v+A0*v+exp(-l)*(A1*v);
nep.resnorm=@(l,v) norm(nep.M(l,v))/((abs(l)+norm(A0,inf)+abs(exp(-l))*norm(A1,inf))*norm(v));




% as error measure the distance between the closest eigenvalue
M_err = @(l,v) min(abs(l-E));
nep.resnorm=M_err;
% solve the projected problem with TIAR
nep.Md =@(j,v) nep.Md(j)*v; % adapt the interface for TIAR
opts.maxit=m; opts.tol=1e-5; opts.disp=1;  k=Inf;
opts.sigma=0;    opts.gamma=1;  opts.p=1;
[ W, E ] = tiar( nep, k, opts );

% plot the eigenvalues approximation extracted by the projected problem
figure(2)
plot(real(E),imag(E),'d'); hold on
theta=linspace(0,2*pi,100);
plot(cos(theta),sin(theta),'k-')

% plot the decay of the q'*B*q coefficients (that are not used practically!)
%figure
%d=diag(Omega);
%semilogy(abs(d))
clc
fprintf("INFLAN %f seconds\n",t1)
fprintf("TIAR %f seconds\n",t2)