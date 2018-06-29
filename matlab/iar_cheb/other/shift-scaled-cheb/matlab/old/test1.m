close all
clear all
clc

%coeff=[4; 2; -1; 5; -2];
coeff=rand(7,1);
n=length(coeff)-1;

p=@(x) polyval(flip(coeff),x);

% Chebyshev polynomials of the first kind
T=@(n,x) cos(n*acos(x));

rho=2;    gamma=1;
% Chebyshev polynomials (shifted and rescaled) of the first kind
TT=@(n,x) T(n,rho*x+gamma);



cc=mon2cheb(rho,gamma,coeff);
c=naive_mon2cheb(TT,coeff);
fprintf("Error=%e\n",norm(c-cc));

xp=linspace(-1,1,100);
yp=zeros(size(xp));
for j=1:length(xp)
    yp(j)=pT(xp(j),c,TT);
end
plot(xp,yp,'--r')


% plot the polynomial in both basis (should be the same)
hold on
for j=1:length(xp)
    yp(j)=p(xp(j));
end
plot(xp,yp,'-k')


function y=pT(x,c,TT)
    y=0;
    for j=0:length(c)-1
        y=y+c(j+1)*TT(j,x);
    end
end



function c=naive_mon2cheb(TT,coeff)
% nodes
p=@(x) polyval(flip(coeff),x);
n=length(coeff)-1;
xx=rand(n+1,1);
A=zeros(n+1,n+1);
b=zeros(n+1,1);
for i=1:n+1
    b(i)=p(xx(i));
    for j=0:n
        A(i,j+1)=TT(j,xx(i));
    end
end
c=A\b;

end

function c=mon2cheb(rho,gamma,a)
    n=length(a)-1;

    % constants in the three term recurrence
    alpha=1/(2*rho);    beta=-gamma/rho;

    b=zeros(n+3,1);
    bb=zeros(n+3,1);
    for j=n:-1:0
        bb(1)=alpha*b(2)+beta*b(1)+a(j+1);
        bb(2)=beta*b(2)+alpha*b(3)+2*alpha*b(1);
        for k=3:n-j-1
            bb(k)=alpha*b(k-1)+beta*b(k)+alpha*b(k+1);
        end
        if n-j>2
            bb(n-j)=alpha*b(n-j-1)+beta*b(n-j);
        end
        if n-j+1>2
            bb(n-j+1)=alpha*b(n-j);
        end
        b=bb;
        bb=0*bb;
    end
    c=b(1:n+1);
end



