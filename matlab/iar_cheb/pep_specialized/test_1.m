close all
clear all
clc


a=-1;   b=2;


n=100; m=10;  d=100;
for j=1:d+1
    coeff{j}=rand(n)/factorial(j);
end
%coeff{1}=eye(n);
nep.coeff=coeff;
nep.A0inv=inv(coeff{1});

nep.MMeval=@(l)  M_pep_eval(coeff,l);
nep.Mdd=@(j) M_pep_der(coeff,j);

nep.M0solver=@(x) nep.MMeval(0)\x;
nep.err=@(lambda,v) norm(nep.MMeval(lambda)*v);
nep.n=n;
nep.a=a;    nep.b=b;    
%v=zeros(n,1);   v(1)=1;
v=ones(n,1);    v=v/norm(v);
%v=rand(n,1);   v=v/norm(v);


[ V, H ] = InfArn_change_basis( nep, v, m ); 
%[ V, H ] = InfArn_dep( nep, v, m ); 

V=V(1:n,:);
[ err, conv_eig_IAR_Chebyshev ] = iar_error_hist( nep, V, H, '--r' );


%figure
%AA=plot(real(conv_eig_IAR),imag(conv_eig_IAR),'d');    hold on
%BB=plot(real(conv_eig_IAR_Chebyshev),imag(conv_eig_IAR_Chebyshev),'*');
%legend([AA, BB],'Taylor','Chebyshev')
%title('Converged Ritz values')



function MM=M_pep_eval(coeff,l)
    MM=zeros(size(coeff{1}));
    for j=1:length(coeff)
        MM=MM+l^(j-1)*coeff{j};
    end
end

function MM=M_pep_der(coeff,j)
    d=length(coeff);
    n=length(coeff{1});
    if j<d
        MM=factorial(j)*coeff{j+1};
    else
        MM=zeros(n);
    end
end