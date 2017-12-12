% converts cheb to monomial

close all
clear all
clc

n=10;
b=rand(n,1);
%b=[5; 3; 2; 6; 4];

aa1=g_cheb2mon(b);
aa2=cheb2mon(b);

norm(aa1-aa2)

%[aa1 aa2]



function a=g_cheb2mon(b)
% a is a column vector

n=length(b);
b=flip(b);
b(end+2)=0;
bb=zeros(n+2,1);
a=zeros(n,1);


for j=1:n 
   

   
   for k=n-j+1:-1:2
       bb(k)=b(k+1)-bb(k+2);
   end
   
   bb(1)=(b(2)-bb(3))/2;
   a(j)=(2^(j-1))*(b(1)-bb(2));
   
   
   b=bb;
   bb=0*bb;
end

a=flip(a);

end


function a = mon2cheb(b)
%CHEB2MON  Monomial to Chebyshev basis conversion.
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
%
%   See also   MON2CHEB, CHEBPOLY

%   Zoltn Csti
%   2015/03/31


% Construct the Chebyshev polynomials of the first kind
N = length(b)-1;
C = chebpoly(N);
% Create the transformation matrix
A = zeros(N+1);
for k = 1:N+1
   A(k:N+1,k) = C{N+2-k};
end
% Perform the basis conversion by solving the linear system
a = A\b(:);
end

function b = cheb2mon(a)
%CHEB2MON  Chebyshev to monomial basis conversion.
%   B = CHEB2MON(A) converts polynomial A given in Chebyshev basis to 
%   monomial basis B. The polynomial must be given with its coefficients in
%   descending order, i.e. A = A_N*T_N(x) + ... + A_1*T_1(x) + A_0*T_0(x)
%
%   Example: 
%    Suppose we have a polynomial in Chebyshev basis: 
%    a2*T_2(x) + a1*T_1(x) + a0*T_0(x), where T_0=1, T_1=x, T_2=2x^2-1
%    and for example a2=1, a1=0, a0=-1.
%    We want to express the polynomial in the monomial base {1,x,x^2), i.e.
%    a2*T_2(x) + a1*T_1(x) + a0*T_0(x) = b2*x^2 + b1*x + b0,
%    where b = [b2 b1 b0] is sought.
%    Solution:
%      a = [1 0 -1];
%      b = cheb2mon(a);
%
%   See also   CHEB2MON, CHEBPOLY

%   Zoltn Csti
%   2015/03/31


% Construct the Chebyshev polynomials of the first kind
N = length(a)-1;
C = chebpoly(N);
% Create the transformation matrix
A = zeros(N+1);
for k = 1:N+1
   A(k:N+1,k) = C{N+2-k};
end
% Perform the basis conversion using a matrix-vector product
b = A*a(:);
end

function C = chebpoly(n)
% C = CHEBPOLY(N) returns Chebyshev polynomials of the first kind from 
% degree 0 to n

C = cell(n+1,1);   % using cell array because of the different matrix sizes
C{1} = 1;          % T_0 = 1
C{2} = [1 0];      % T_1 = x
for k = 3:n+1      % T_n = 2xT_n-1 - T_n-2
    C{k} = [2*C{k-1} 0] - [0 0 C{k-2}];
end
end
