close all
clear all
clc

maxit=100;
n=1000; 
m=2*maxit+2;

% matrices defining the problem
A{1}=spdiags(rand(n,3),-1:1,n,n)+spdiags(rand(n,3),-1:1,n,n)*1i;        
A{1}=A{1}+A{1}';

A{2}=spdiags(rand(n,3),-1:1,n,n)+spdiags(rand(n,3),-1:1,n,n)*1i;        
A{2}=A{2}+A{2}';

A{3}=spdiags(rand(n,3),-1:1,n,n)+spdiags(rand(n,3),-1:1,n,n)*1i;        
A{3}=A{3}+A{3}';

A{4}=spdiags(rand(n,3),-1:1,n,n)+spdiags(rand(n,3),-1:1,n,n)*1i;        
A{4}=A{4}+A{4}';
A{4}=1e-12*A{4};

A{5}=spdiags(rand(n,3),-1:1,n,n)+spdiags(rand(n,3),-1:1,n,n)*1i;        
A{5}=A{5}+A{5}';

% functions defining the problem
f{1}=@(l) -l; 
f{2}=@(l) eye(size(l));
f{3}=@(l) l;
%f{4}=@(l) sqrtm(10*l+100*eye(size(l)));
f{4}=@(l) l^3;
f{5}=@(l) expm(-l);

mu=0;
nep=setup_nep(f,A,mu,m,[]);



Q1=rand(n,1)+rand(n,1)*1i; Q1=Q1/norm(Q1);
[Q,T,omega] = Lanczos_nep_full_complex(nep,maxit-1,Q1);
[Q,~]=qr(Q,0);

% set up the projected nep
proj_A=cell(1,length(A));
for j=1:length(A)
    proj_A{j}=Q'*(A{j}*Q);
end
proj_nep=setup_nep(f,proj_A,mu,m,[]);

% solve the original problem with TIAR
opts.maxit=maxit; opts.tol=1e-10; opts.disp=1; opts.v0=rand(n,1); opts.p=Inf; 
[ V, E ] = tiar( nep, 200, opts );

% solve the projected problem with TIAR
clear opts
%proj_nep.resnorm=@(l,v) min(abs(l-E));
proj_nep.resnorm=@(l,v) nep.resnorm(l,Q*v);
opts.maxit=maxit; opts.tol=1e-6; opts.disp=1; opts.p=1;
[ ~, EE ] = tiar( proj_nep, 100, opts );

figure(5)
AA=plot(real(EE),imag(EE),'ko'); hold on
BB=plot(real(E),imag(E),'b*');
legend('Lanczos','TIAR')
