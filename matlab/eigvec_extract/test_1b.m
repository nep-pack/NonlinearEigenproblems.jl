%% TEST ON A DEP
close all
clear all
clc

%% generate the NEP
rng('default')
rng(1)
n=100; A0 = spdiags(rand(n,3), -1:1, n, n);
A1 = spdiags(rand(n,3), -1:1, n, n);
nep.M =@(l,v)  -l*v+A0*v+exp(-l)*(A1*v);

AA=cell(3,1);
AA{1}=@(v) v+A1*v;
AA{2}=@(v) A1*v;
AA{3}=@(v) A1*v;
nep.Md=@(j,v) AA{min(j,3)}(v)*(-1)^j;
[L,U,P,Q] = lu(A0+A1);

nep.M_inv=@(v) Q*(U\(L\(P*v)));
nep.n=n;
nep.resnorm=@(l,v) norm(nep.M(l,v))/((abs(l)+norm(A0,inf)+abs(exp(-l))*norm(A1,inf))*norm(v));
opts.maxit=50; opts.tol=1e-10; opts.disp=0; opts.v0=ones(n,1); opts.p=1; k=inf;


%% solve the NEP with TIAR
[ W, E ] = tiar( nep, k, opts ); 
vv=W(:,1);  

% check we have an eigenpair
%nep.resnorm(E(1),vv)
%norm(nep.M(E(1),vv))

%% extraction of the eigenvector/nullspace of A
% we want to compute the kernel of this matrix
A=-E(1)*speye(n)+A0+exp(-E(1))*A1;
norm(A*vv)


x0=rand(n,1);
%x0=vv;
options = optimoptions(@fmincon,'Algorithm','interior-point','MaxIterations',1500,'MaxFunctionEvaluations',1e6,'Display','iter','StepTolerance',1e-14);
[x,fval,exitflag,output]=fmincon(@(x) norm(A*x)^2, x0, [], [], ones(1,n),1,[],[],[],options);
norm(A*x)
nep.resnorm(E(1),x)
