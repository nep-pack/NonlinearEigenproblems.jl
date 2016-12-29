% TESTS ON THE BASIS MATRIX AFTER RESTART 
clear all
close all
clc

% TEST PROBLEM
n=3;
I=speye(n);
e = ones(n,1);  
A2 = spdiags([e -2*e e], -1:1, n, n);     
A1 = spdiags([e e e], -1:1, n, n); 
%A1=rand(n,n);   A2=rand(n,n);   A1=sparse(A1);  A2=sparse(A2);

n=length(A1);
nep.n=n;

% DEFINITION OF THE PROBLEM
%nep.M =@(lambda) -lambda*I + A2 + exp(-lambda)*A1;

[L,U,P,Q,R] = lu(A1+A2);    R=diag(R);

nep.M0solver=@(v) Q*(U\(L\(P*(v./R))));
nep.err=@(lambda,v)  norm(-lambda*v + A2*v + exp(-lambda)*A1*v)/...
                    (abs(lambda)+ norm(A2,1)+abs(exp(-lambda))*norm(A1,1))*norm(v);
    
nep.Md=@(j,v)   (j==0)*(A1*v+A2*v)+...
                (j==1)*(-v-A1*v)+...
                (j>1)*((-1)^j)*A1*v;
                
            
            

v=ones(n,1);
v=v/norm(v);

m=20;
E=zeros(m,1);
for i=1:m
    v=B(nep,v);
    lambda=norm(v);
    v=v/lambda;
    E(i)=nep.err(1/lambda,v(1:n,:)); 
    lambda    
end
   



semilogy(E,'--*')