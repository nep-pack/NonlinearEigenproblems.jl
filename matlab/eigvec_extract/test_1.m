%% TEST ON A DEP
close all
clear all
clc

%% generate the NEP
rng('default')
rng(1)
n=5000; A0 = spdiags(rand(n,3), -1:1, n, n);
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


% starting vector
%b=randn(n,1);

% choose this option to start with an eigenvector approximation
b=vv+randn(n,1)*1e-6;

nv1=[]; nv2=[]; nv3=[];

% generate an orthogonal basis of the Krylov space
m=200;  [ V, H ] = arnoldi( A, b, m+1 );

% Solve the problem with Arnoldi
for k=1:m    
    VV=V(:,1:k);    HH=H(1:k,1:k);    
    [Z,D] = eig(HH);    D=diag(D);     
    [~,idx]=sort(abs(D));  v=VV*Z(:,idx(1));        v=v/norm(v);
    nv1=[nv1,nep.resnorm(E(1),v)];
end
fprintf("Residual with Arnoldi %e\n",norm(A*v))

% Solve the problem with Harmonic Arnoldi
for k=1:m    
    VV=V(:,1:k);    HH=H(1:k+1,1:k);    
    [Z,D]=eig(H(1:k,1:k),HH'*HH); D=diag(D);  D=1./D;
    [~,idx]=sort(abs(D));  v=VV*Z(:,idx(1));        v=v/norm(v);
    nv2=[nv2,nep.resnorm(E(1),v)];
end
fprintf("Residual with Harmonic Arnoldi %e\n",norm(A*v))

% Solve the problem with GMRES
for k=1:m    
    VV=V(:,1:k);        HH=H(1:k+1,1:k);    
    z=HH\eye(k+1,1);    v=VV*z;        v=v/norm(v);
    nv3=[nv3,nep.resnorm(E(1),v)];
end
fprintf("Residual with gmres %e\n",norm(A*v))

% plot results
close all
AA=semilogy(nv1,'b');  hold on
BB=semilogy(nv2,'k');
CC=semilogy(nv3,'g');
legend([AA BB CC],'Arnoldi','Harmonic Arnoldi','gmres');