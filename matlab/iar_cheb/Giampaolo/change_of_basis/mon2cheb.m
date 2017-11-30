function bb=mon2cheb(a)
%mon2cheb  Monomial to Chebyshev basis conversion.
%   A = MON2CHEB(B) converts polynomial B given in monomial basis to 
%   Chebyshev basis A. The polynomial must be given with its coefficients
%   in descending order, i.e. B = B_N*x^N + ... + B_1*x + B_0
%
%   Example: 
%    Suppose we have a polynomial in the monomial basis: 
%    b2*x^2 + b1*x + b0, 
%    with b2=2, b1=0, b0=-2 for example.
%    We want to express the polynomial in the Chebyshev base 
%    {T_0(x),T_1(x),T_2(x)}, where T_0=1, T_1=x, T_2=2x^2-1, i.e.
%    a2*T_2(x) + a1*T_1(x) + a0*T_0(x) = b2*x^2 + b1*x + b0,
%    where a = [a2 a1 a0] is sought.
%    Solution:
%      b = [2 0 -2];
%      a = mon2cheb(b);

n=length(a);
a=flipud(a);
bb=zeros(n+2,1);
b=0*bb;


for j=n:-1:1 
    b(1)=a(j)+bb(2)/2;
    b(2)=bb(1)+bb(3)/2;
    for k=3:n-j+2
        b(k)=(bb(k-1)+bb(k+1))/2;
    end
    k=n-j+3;
	b(k)=b(k-1)/2;
    
    bb=b;   
end

bb=bb(1:n);
bb=flipud(bb);
end