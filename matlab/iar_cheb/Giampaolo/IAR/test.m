close all
clear all
clc

n=4;
A0=rand(n); A1=rand(n); I=eye(n);

nep.MMeval=@(l)  -l^2*I + A0 + A1*exp(-l);
nep.Mdd=@(j)                           ...
                (j==0)*(A0 + A1) +    ...
                (j==1)*(-A1) +          ...
                (j==2)*(-2*I+A1) +      ...
                (j>2)*((-1)^j*A1);
nep.M0solver=@(x) nep.MMeval(0)\x;
nep.err=@(lambda,v) norm(nep.MMeval(lambda)*v);
nep.n=n;

v=zeros(n,1);   v(1)=1;
m=20;
[ V, H ] = InfArn( nep, v, m ); V=V(1:n,:);
[ err, conv_eig_IAR ] = iar_error_hist( nep, V, H, '-k' );