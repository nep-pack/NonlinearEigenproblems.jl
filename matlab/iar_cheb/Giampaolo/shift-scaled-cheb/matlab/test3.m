close all
clear all
clc

n=10;
%c=1./exp(1:1:n);
%c=n:-1:1;
c=rand(n,1);

a=-1;   b=0;
rho=2/(b-a); gamma=(a+b)/(a-b);

a=cheb2mon(rho,gamma,c);
cc=mon2cheb(rho,gamma,a);

double(norm(c-cc))

C=semilogy(abs(c),'--k');
hold on
A=semilogy(abs(a),'-r');

legend([C, A],"Chebishev","Monomials")