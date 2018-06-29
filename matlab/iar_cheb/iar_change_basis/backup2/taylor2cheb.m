function c=taylor2cheb(rho,gamma,a)
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
%        bb(1)=alpha*b(2)+beta*b(1)+a(j+1);
        bb(1)=alpha*b(2)+beta*b(1)+a(j+1)*factorial(j+1);
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
    c=c(:,1);
end