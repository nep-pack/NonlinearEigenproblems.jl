close all
clear all
clc


a=-4;   b=2;


n=100; 
I=eye(n); A0=rand(n); A1=rand(n); m=100; tau1=2;
A2=rand(n); tau2=1.4;

nep.MMeval=@(l)  -l*I + A0 + A1*exp(-l*tau1) + A2*exp(-l*tau2);
nep.Mdd=@(j)                                    ...
                (j==0)*(A0+A1+A2) +             ...
                (j==1)*(-I-tau1*A1-tau2*A2) +   ...
                (j>1)*((-tau1)^j*A1+(-tau2)^j*A2);

nep.M0solver=@(x) nep.MMeval(0)\x;
nep.err=@(lambda,v) norm(nep.MMeval(lambda)*v);
nep.n=n;
nep.A0=A0;  nep.A1=A1;  nep.A2=A2;  nep.I=I;
nep.a=a;    nep.b=b;    nep.tau1=tau1;  nep.tau2=tau2;

%v=zeros(n,1);   v(1)=1;
v=ones(n,1);    v=v/norm(v);
%v=rand(n,1);


[ V, H ] = InfArn_change_basis( nep, v, m ); 
%[ V, H ] = InfArn_dep( nep, v, m ); 

V=V(1:n,:);
[ err, conv_eig_IAR_Chebyshev ] = iar_error_hist( nep, V, H, '--r' );


%figure
%AA=plot(real(conv_eig_IAR),imag(conv_eig_IAR),'d');    hold on
%BB=plot(real(conv_eig_IAR_Chebyshev),imag(conv_eig_IAR_Chebyshev),'*');
%legend([AA, BB],'Taylor','Chebyshev')
%title('Converged Ritz values')

