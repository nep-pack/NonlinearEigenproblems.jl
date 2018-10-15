close all
clear all
clc

maxit=5;
n=3; 
m=2*maxit+2;

% matrices defining the problem
A{1}=[ 1 -4 -2
    -4  6  3
    -2  3  7];
A{2}=[ 3 5 4
     5 1 -1
     4 -1 11];
A{3}=[ 2   1 -11
     1   2  3
    -11 3 -3];
A{4}=[ -1  1   12
      1  -4  -2
     12  -2  7];


% functions defining the problem
f{1}=@(l) eye(size(l));
f{2}=@(l) -l;
f{3}=@(l) expm(-l);
f{4}=@(l) sqrtm(eye(size(l))-2*l);


mu=0;
nep=setup_nep(f,A,mu,m,[]);



Q1=ones(n,1); Q1=Q1/norm(Q1);
[Q,T,omega] = Lanczos_nep_full(nep,maxit-1,Q1);
Q


return
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
