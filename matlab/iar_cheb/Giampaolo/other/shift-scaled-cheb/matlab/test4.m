close all
clear all
clc

% set scaling and shift
%rho=2;    gamma=3;
%rho=1;    gamma=0;
%a=-3;   b=3;
a=-1;   b=1;
rho=2/(b-a); gamma=(a+b)/(a-b);

n=50;
I=eye(n);
P=zeros(n); Pinv=zeros(n);
for j=1:n
    P(:,j)=cheb2mon(rho,gamma,I(:,j));
    Pinv(:,j)=mon2cheb(rho,gamma,I(:,j));
end


c=rand(n,1);
a=P*c;
cc=Pinv*a;
%cc=(Pinv*P)*c;

norm(c-cc)



for j=1:n
    D(j)=1/factorial(j);
    Dinv(j)=factorial(j);
end
 
D=diag(D);
Dinv=diag(Dinv);

PP=D*P;
aa=PP*c;
a=Dinv*aa;

PPinv=Dinv*Pinv;
cc=PPinv*a;
cc=D*cc;
norm(c-cc)


% norm(P)
% norm(PP)


%[a,FLAG,RELRES,ITER,resvec] = gmres(Pinv,c,1,1e-12,100,P);
%semilogy(resvec)
%norm(c-cc)


