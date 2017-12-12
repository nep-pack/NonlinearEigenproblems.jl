close all
clear all
clc

% set scaling and shift
%rho=2;    gamma=3;
%rho=1;    gamma=0;
%a=-3;   b=3;
a=-1;   b=0;
rho=2/(b-a); gamma=(a+b)/(a-b);

n=10;
I=eye(n);
P=zeros(n); Pinv=zeros(n);
for j=1:n
    P(:,j)=cheb2mon(rho,gamma,I(:,j));
    Pinv(:,j)=mon2cheb(rho,gamma,I(:,j));
end

P
%Pinv

%P*Pinv

LL = @(k) L(k)*(nep.b-nep.a)/4;
A=L(n)';
A

%K=krylov(L(n),I(:,1))



function  K=krylov(A,v)
    n=length(v);
    K=zeros(n);
    K(:,1)=v;
    for j=1:n-1
        K(:,j+1)=A*K(:,j);
    end
end


function Lmat=L(k)
    if k==1
        Lmat=2;
    else
        Lmat=diag([2, 1./(2:k)])+diag(-1./(1:(k-2)),-2);
    end
end