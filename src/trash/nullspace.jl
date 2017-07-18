workspace()

A=sparse([1 1.0 1.0 ; 0 1 1; -1 0 0]);
F=lufact(A)
F[:L]*F[:U] - (F[:Rs] .* A)[F[:p], F[:q]]
UU=F[:U]; UU[end,end]=1;
n=size(A,1);

b=zeros(n);
b[end]=1;
x=UU\b;

norm(A*x)
