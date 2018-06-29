close all
clear all
clc

T=@(i,x) cos(i*acos(x));
a=-1;   b=1;
c=(a+b)/(a-b);
k=2/(b-a);

N=10;
n=20;
tau=1;
k=2/(b-a);



y=rand(n,N+1);
x=rand(n,N);

y0=zeros(n,1);
for i=1:N
    y0=y0+T(i,c)*y(:,i+1);
end
y0=-y0;

for i=0:N-1
    y0=y0+T(i,-k*tau+c)*x(:,i+1);
end