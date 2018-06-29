close all
clear all
clc


a=-1;   b=0;

n=1000; 
m=100;



I=eye(n);


D=10*rand(n,1);
[Q,~]=qr(rand(n));
A=Q*diag(D)*Q';

nep.MMeval=@(l)  A-l*I;
nep.Mdd=@(j)    (j==0)*(A)+(j==1)*(-I);

            

nep.M0solver=@(x) nep.MMeval(0)\x;
nep.err=@(lambda,v) norm(nep.MMeval(lambda)*v);
nep.n=n;
nep.a=a;    nep.b=b;

%v=zeros(n,1);   v(1)=1;
%v=ones(n,1);    v=v/norm(v);
v=rand(n,1);


[ V, H ] = InfArn_change_basis( nep, v, m ); 
V=V(1:n,:);
[ err, conv_eig_IAR_Chebyshev ] = iar_error_hist( nep, V, H, '--r' );


%figure
%AA=plot(real(conv_eig_IAR),imag(conv_eig_IAR),'d');    hold on
%BB=plot(real(conv_eig_IAR_Chebyshev),imag(conv_eig_IAR_Chebyshev),'*');
%legend([AA, BB],'Taylor','Chebyshev')
%title('Converged Ritz values')

