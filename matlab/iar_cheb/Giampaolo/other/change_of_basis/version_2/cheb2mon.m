function a=cheb2mon(b)
%cheb2mon:  Monomial to Chebyshev basis conversion.
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



n=length(b);
b=flipud(b);
b(end+2)=0;
bb=zeros(n+2,1);
a=zeros(n,1);


for j=1:n 
   

   
   for k=n-j+1:-1:2
       bb(k)=2*b(k+1)-bb(k+2);
   end
   
   bb(1)=b(2)-bb(3)/2;
   a(j)=b(1)-bb(2)/2;
   
   
   b=bb;
   bb=0*bb;
end

a=flipud(a);

end
