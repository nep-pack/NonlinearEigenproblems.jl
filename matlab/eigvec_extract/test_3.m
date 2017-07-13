%% TEST ON A RANDOM SINGULAR MATRIX
close all
clear all
clc

%% generate the problem
n=1000;
A=randn(n);  

d=eig(A);
A=A-eye(length(A))*(d(1)+eps*100);

%% extraction of the eigenvector/nullspace of A
% we want to compute the kernel of this matrix

% starting vector
b=randn(n,1);

% choose this option to start with an eigenvector approximation
%b=vv+randn(n,1)*1e-6;

nv1=[]; nv2=[]; nv3=[];

% generate an orthogonal basis of the Krylov space
m=200;  [ V, H ] = arnoldi( A, b, m+1 );

% Solve the problem with Arnoldi
for k=1:m    
    VV=V(:,1:k);    HH=H(1:k,1:k);    
    [Z,D] = eig(HH);    D=diag(D);     
    [~,idx]=sort(abs(D));  v=VV*Z(:,idx(1));        v=v/norm(v);
    nv1=[nv1, norm(A*v)];
end
fprintf("Residual with Arnoldi %e\n",norm(A*v))

% Solve the problem with Harmonic Arnoldi
for k=1:m    
    VV=V(:,1:k);    HH=H(1:k+1,1:k);    
    [Z,D]=eig(H(1:k,1:k),HH'*HH); D=diag(D);  D=1./D;
    [~,idx]=sort(abs(D));  v=VV*Z(:,idx(1));        v=v/norm(v);
    nv2=[nv2,norm(A*v)];
end
fprintf("Residual with Harmonic Arnoldi %e\n",norm(A*v))

% Solve the problem with GMRES
for k=1:m    
    VV=V(:,1:k);        HH=H(1:k+1,1:k);    
    z=HH\eye(k+1,1);    v=VV*z;        v=v/norm(v);
    nv3=[nv3,norm(A*v)];
end
fprintf("Residual with gmres %e\n",norm(A*v))


% plot results
close all
AA=semilogy(nv1,'b');  hold on
BB=semilogy(nv2,'k');
CC=semilogy(nv3,'g');
legend([AA BB CC],'Arnoldi','Harmonic Arnoldi','gmres');

%% result of one iteration of residual inverse iteration
v=A\b;  v=v/norm(v);
fprintf("One iteration of inverse iteration gives error %e\n",norm(A*v))