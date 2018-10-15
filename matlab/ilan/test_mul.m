close all
clear all
clc
n=10;
k=4;

Z=rand(n,k)+rand(n,k)*1i;
W=rand(n,k)+rand(n,k)*1i;

gamma=sum(sum(conj(Z).*W));%=Z(:)'*W(:);
gamma2=(Z(:)')*W(:);
norm(gamma-gamma2)