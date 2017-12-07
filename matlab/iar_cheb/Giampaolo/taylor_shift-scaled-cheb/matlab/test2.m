close all
clear all
clc

n=50;
c=1./exp(1:1:n);
c=rand(n,1);

rho=1;  gamma=0;
a=cheb2mon(rho,gamma,c);
cc=mon2cheb(rho,gamma,a);
norm(c-cc)

C=semilogy(c,'--k');
hold on
A=semilogy(a,'-r');

legend([C, A],"Chebishev","Monomials")