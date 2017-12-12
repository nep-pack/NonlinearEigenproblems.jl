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

p=@(x) Taylor_polyval(a,x);
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
c=c(:,1);

end



function y=Taylor_polyval(a,x)
    y=0;
    for j=1:length(a)
        y=y+(a(j)/factorial(j))*x^(j-1);
    end
end