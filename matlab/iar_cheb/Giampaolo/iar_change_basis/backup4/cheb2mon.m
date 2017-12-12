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
    a=zeros(n+3,1);     b=zeros(n+3,1);    
    bb=zeros(n+3,1);    bb(1:n+1)=c;   
    %bb=sym(bb); b=sym(b);   a=sym(a);    
    for j=1:1:n+1
        for k=n-j+1:-1:2
            b(k)=(bb(k+1)-beta*b(k+1)-alpha*b(k+2))/alpha;
        end
        b(1)=(bb(2)-beta*b(2)-alpha*b(3))/(2*alpha);
        a(j)=bb(1)-alpha*b(2)-beta*b(1);
        
        bb=b; b=0*b;
    end
    a=a(1:n+1,1);
end