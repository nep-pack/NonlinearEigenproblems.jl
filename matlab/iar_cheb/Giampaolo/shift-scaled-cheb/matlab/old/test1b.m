close all
clear all
clc

%coeff=[4; 2; -1; 5; -2];
a=rand(6,1);
n=length(a)-1;

p=@(x) polyval(flip(a),x);

% Chebyshev polynomials of the first kind
T=@(n,x) cos(n*acos(x));

rho=2;    gamma=1;
% Chebyshev polynomials (shifted and rescaled) of the first kind
TT=@(n,x) T(n,rho*x+gamma);



cc=mon2cheb(rho,gamma,a);
c=naive_mon2cheb(rho,gamma,a);
aa=cheb2mon(rho,gamma,c);
[a aa]
%coeff-a
fprintf("Error=%e\n",norm(c-cc));

a=cheb2mon(2,3,[2 -4 -1])

x=linspace(-1,1,100);
y=cheb_polyval(x,c,rho,gamma);

hold on
y=polyval(flip(a),x);
plot(x,y,'-k')


function y=cheb_polyval(x,c,rho,gamma)
% cheb_polyval: evaluate a polynomial written in 
% shifted-and-scaled Chebyshev basis
% Given the polynomial
% p(x)=\sum_{j=0}^n c(j+1) T_j(rho*x+gamma)
% this function returns the vector y such that y(i)=p(x(i)) for every 
% component of x

    % Chebyshev polynomials of the first kind
    T=@(n,x) cos(n*acos(x));

    % Chebyshev polynomials (shifted and rescaled) of the first kind
    TT=@(n,x) T(n,rho*x+gamma);
    y=zeros(size(x));
    for i=1:length(x)
        y(i)=0;
        for j=0:length(c)-1
            y(i)=y(i)+c(j+1)*TT(j,x(i));
        end
    end
end



function c=naive_mon2cheb(rho,gamma,a)
% DO NOT USE THIS FUNCTION! ONLY FOR VALIDATION!
% THIS FUNCTION IT IS BASED ON INTERPOLATION AND IT IS EXPECTED TO NOT BE
% NUMERICALLY STABLE FOR COEFF THAT IS A LONG VECTOR.
%naive_cheb2mon: shifted-and-scaled Chebyshev basis conversion to Monomial basis.
%   a = cheb2mon(rho,gamma,c) converts a polynomial written in 
%   shifted-and-scaled Chebyshev basis
%   p=\sum_{j=0}^n c(j+1) T_j(rho*x+gamma)
%   to a polynomial written in monomial basis
%   p=\sum_{j=0}^n a(j+1) x^j
%   where T_j(x) is the j-th Chebyshev polynomial of the first kind. 
%
%   Example: 
%    Suppose we have a polynomial in the shifted-and-scaled Chebyshev basis: 
%    2*T_0(2x+3)-4*T_1(2x+3)-T_2(2x+3) 
%    This polynomial is expressed in monomial basis as
%    a(1)+a(2)*x+a(3)*x^2
%    where 
%    a=cheb2mon(2,3,[2 -4 -1])

% Chebyshev polynomials of the first kind
T=@(n,x) cos(n*acos(x));

% Chebyshev polynomials (shifted and rescaled) of the first kind
TT=@(n,x) T(n,rho*x+gamma);

p=@(x) polyval(flip(a),x);
n=length(a)-1;
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
%cheb2mon: Monomial basis to shifted-and-scaled Chebyshev basis conversion.
%   c = cheb2mon(rho,gamma,a) converts a polynomial written in 
%   monomial basis
%   p(x)=\sum_{j=0}^n a(j+1) x^j
%   to a polynomial written in shifted-and-scaled Chebyshev basis
%   p(x)=\sum_{j=0}^n c(j+1) T_j(rho*x+gamma)
%   where T_j(x) is the j-th Chebyshev polynomial of the first kind. 
%
%   Example: 
%    Suppose we have a polynomial in the monomial basis: 
%    2-4*x-x^2 
%    This polynomial is expressed in shifted-and-scaled Chebyshev basis as:
%    c(1)*T_0(2*x+3)+c(2)*T_1(2*x+3)+a(3)*T_2(2*x+3)
%    where 
%    c=cheb2mon(2,3,[2 -4 -1])

    n=length(a)-1;

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




function a=cheb2mon(rho,gamma,c)
%cheb2mon: shifted-and-scaled Chebyshev basis conversion to Monomial basis.
%   a = cheb2mon(rho,gamma,c) converts a polynomial written in 
%   shifted-and-scaled Chebyshev basis
%   p(x)=\sum_{j=0}^n c(j+1) T_j(rho*x+gamma)
%   to a polynomial written in monomial basis
%   p(x)=\sum_{j=0}^n a(j+1) x^j
%   where T_j(x) is the j-th Chebyshev polynomial of the first kind. 
%
%   Example: 
%    Suppose we have a polynomial in the shifted-and-scaled Chebyshev basis: 
%    2*T_0(2x+3)-4*T_1(2x+3)-T_2(2x+3) 
%    This polynomial is expressed in monomial basis as
%    a(1)+a(2)*x+a(3)*x^2
%    where 
%    a=cheb2mon(2,3,[2 -4 -1])

    n=length(c)-1;
    alpha=1/(2*rho);    beta=-gamma/rho;
    a=zeros(n+1,1);     b=zeros(n+2,1);    
    bb=zeros(n+2,1);    bb(1:n+1)=c;   
    for j=1:1:n+1
        for k=n-j+1:-1:2
            b(k)=(bb(k+1)-beta*b(k+1)-alpha*b(k+2))/alpha;
        end
        b(1)=(bb(2)-beta*b(2)-alpha*b(3))/(2*alpha);
        a(j)=bb(1)-alpha*b(2)-beta*b(1);
        bb=b; b=0*b;
    end
end