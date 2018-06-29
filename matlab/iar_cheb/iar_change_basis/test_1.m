close all
clear all
clc


a=-1;   b=1;


n=100; A0=rand(n); A1=rand(n); m=5;

n=4;

A0=[1 2 5 9
    4 8 3 8
    1 5 -4 5
    1 3 5 4
    ];


A1=[4 2 5 9
    5 -8 3 8
    1 -5 -3 5
    4 3 5 4
    ];

m=5;

I=eye(n);


nep.MMeval=@(l)  -l*I + A0 + A1*exp(-l);
nep.Mdd=@(j)                            ...
                (j==0)*(A0 + A1) +      ...
                (j==1)*(-I-A1) +        ...
                (j>1)*((-1)^j*A1);

            

nep.M0solver=@(x) nep.MMeval(0)\x;
nep.err=@(lambda,v) norm(nep.MMeval(lambda)*v);
nep.n=n;
nep.A0=A0;  nep.A1=A1;  nep.I=I;
nep.a=a;    nep.b=b;

%v=zeros(n,1);   v(1)=1;
v=ones(n,1);    v=v/norm(v);
%v=rand(n,1);


[ V, H ] = InfArn_change_basis( nep, v, m ); 
V=V(1:n,:);
[ err, conv_eig_IAR_Chebyshev ] = iar_error_hist( nep, V, H, '--r' );


%figure
%AA=plot(real(conv_eig_IAR),imag(conv_eig_IAR),'d');    hold on
%BB=plot(real(conv_eig_IAR_Chebyshev),imag(conv_eig_IAR_Chebyshev),'*');
%legend([AA, BB],'Taylor','Chebyshev')
%title('Converged Ritz values')

