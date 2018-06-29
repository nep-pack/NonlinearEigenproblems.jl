close all
clear all
clc

n=4;
a=-1;   b=0;

%A0=rand(n); A1=rand(n); 
A0=[  3  -6   0 4  ;
     -3   4  -8 19 ;
      1 -16 -13 0  ;
    -14  -9   2 9  ]/10;

A1=[  8   2  -13 -3  ;
     -11  9   12  5 ;
      5   2  -16 -13  ;
      7   4   -4   0]/10;

%n=100; A0=rand(n); A1=rand(n); m=100;
I=eye(n);

nep.MMeval=@(l)  -l^2*I + A0 + A1*exp(-l);
nep.Mdd=@(j)                            ...
                (j==0)*(A0 + A1) +      ...
                (j==1)*(-A1) +          ...
                (j==2)*(-2*I+A1) +      ...
                (j>2)*((-1)^j*A1);
nep.M0solver=@(x) nep.MMeval(0)\x;
nep.err=@(lambda,v) norm(nep.MMeval(lambda)*v);
nep.n=n;
nep.A0=A0;  nep.A1=A1;  nep.I=I;
nep.a=a;    nep.b=b;

%v=zeros(n,1);   v(1)=1;
%v=ones(n,1);    v=v/norm(v);
v=rand(n,1);

m=100;

[ V, H ] = InfArn( nep, v, m, 'Taylor' ); 
V=V(1:n,:);
[ err, conv_eig_IAR ] = iar_error_hist( nep, V, H, '-k' );


[ V, H ] = InfArn( nep, v, m, 'Chebyshev' ); 
V=V(1:n,:);
[ err, conv_eig_IAR_Chebyshev ] = iar_error_hist( nep, V, H, '--r' );

% artificially add a legend
AA=semilogy(1e18*rand(2,1),1e18*rand(2,1),'-k');
BB=semilogy(1e18*rand(2,1),1e18*rand(2,1),'--r');
legend([AA, BB],'Taylor','Chebyshev')
title('Convergence history')

figure
AA=plot(real(conv_eig_IAR),imag(conv_eig_IAR),'d');    hold on
BB=plot(real(conv_eig_IAR_Chebyshev),imag(conv_eig_IAR_Chebyshev),'*');
legend([AA, BB],'Taylor','Chebyshev')
title('Converged Ritz values')

