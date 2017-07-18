workspace()

n=20
Q,R=qr(rand(n,n))
D=rand(n)
D[1]=0
A=Q*diagm(D)*Q'
A=sparse(A)

F=lufact(A)
F[:L]*F[:U] - (F[:Rs] .* A)[F[:p], F[:q]]
UU=F[:U]; UU[end,end]=1;
n=size(A,1);

b=zeros(n);
b[end]=1;
x=UU\b;

norm(A*x)
