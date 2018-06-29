%FIXED BY REVERSE ENGINEERING THE ELIAS FORMULA. WRITE IT BETTER.

clc
%n=3;
%%M=randn(n); C=randn(n); K=randn(n);
%A0=randn(n); A1=randn(n);
%
%A0=A0; A1=A1;

n=4;
global A0 A1 a b

a=-1; b=0; 

A0=[0.3000   -0.6000         0    0.4000
   -0.3000    0.4000   -0.8000    1.9000
    0.1000   -1.6000   -1.3000         0
   -1.4000   -0.9000    0.2000    0.9000];

A1=[0.8000    0.2000   -1.3000   -0.3000
   -1.1000    0.9000    1.2000    0.5000
    0.5000    0.2000   -1.6000   -1.3000
    0.7000    0.4000   -0.4000         0];
tau=1;

% Giampaolo edit: change matrices for debug. 
%A0=eye(4);  A1=eye(4);
%A1=0*A1;    % fix the part for A1


%  Compute a very exact solution 
k=2/(b-a) ; c=(a+b)/(a-b);
cheb_vect=@(t,Z) cos((0:(size(Z,2)-1))*acos(t))'; 
cheb2_vect_m1=@(Z)  (0:(size(Z,2)-1))';
Mterm=@(t,X) k*(X*((0:(size(X,2)-1))'.*cheb2_vect_m1(X)));

% QDEP: 
y0comp=@(X,Y) (A0+A1)\(Mterm(c,X)-A0*Y*cheb_vect(c,Y)-A1*Y*cheb_vect(-k*tau+c,Y));


N=3;
x=ones(n,N);
y=ones(n,N+1);

gy0=g_compute_y0(x,y);

ey0=y0comp(x,y);


should_be_equal_columns=[gy0 ey0]
should_be_zero=norm(gy0-ey0)



% function for compuing y0 for this specific DEP
function y0=g_compute_y0(x,y)
    tt=1;
    global A0 A1 a b
    k=2/(b-a) ; c=(a+b)/(a-b);

    
    T=@(n,x) chebyshevT(n,x);   % Chebyshev polynomials of the first kind
    U=@(n,x) chebyshevU(n,x);   % Chebyshev polynomials of the second kind
    
    n=length(A0);
    N=size(x,2);
    
    y0=zeros(n,1);
    for i=1:N-1
        y0=y0+(2*i/a)*U(i-1,1)*x(:,i+1);
    end
    
    for i=1:N+1
        y0=y0+A0*T(i,1)*y(:,i);
    end
    
    for i=1:N+1
        y0=y0+A1*T(i-1,1+2*tt/a)*y(:,i);
    end
    y0=-(A0+A1)\y0;
    
end