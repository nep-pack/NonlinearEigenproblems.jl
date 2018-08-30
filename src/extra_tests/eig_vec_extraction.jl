using LinearAlgebra

n=5;
A=rand(n,n);
A[:,1]=0;

M=zeros(n+1,n);
#M[1:n,:]=A;
#M[n+1,:]=ones(n,1);
A=sparse(A);

M=[A; rand(eltype(A),1,n)]


b=[zeros(eltype(A),n); 1];
#b=A*x;
x=M\b;
y=x[1:n];
println("error=",norm(A*y))
