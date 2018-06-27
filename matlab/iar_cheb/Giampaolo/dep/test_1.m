close all
clear all
clc


a=-1;   b=1;


n=100; 
I=eye(n); A0=rand(n); A1=rand(n); m=100; tau=1;

nep.MMeval=@(l)  -l*I + A0 + A1*exp(-l*tau);
nep.Mdd=@(j)                            ...
                (j==0)*(A0 + A1) +      ...
                (j==1)*(-I-tau*A1) +    ...
                (j>1)*((-tau)^j*A1);

nep.M0solver=@(x) nep.MMeval(0)\x;
nep.err=@(lambda,v) norm(nep.MMeval(lambda)*v);
nep.n=n;
nep.A0=A0;  nep.A1=A1;  nep.I=I;
nep.a=a;    nep.b=b;    nep.tau=tau;

%v=zeros(n,1);   v(1)=1;
v=ones(n,1);    v=v/norm(v);
%v=rand(n,1);


%[ V, H ] = InfArn_change_basis( nep, v, m ); 
[ V, H ] = InfArn_dep( nep, v, m ); 

V=V(1:n,:);
[ err, conv_eig_IAR_Chebyshev ] = iar_error_hist( nep, V, H, '--r' );


%figure
%AA=plot(real(conv_eig_IAR),imag(conv_eig_IAR),'d');    hold on
%BB=plot(real(conv_eig_IAR_Chebyshev),imag(conv_eig_IAR_Chebyshev),'*');
%legend([AA, BB],'Taylor','Chebyshev')
%title('Converged Ritz values')

