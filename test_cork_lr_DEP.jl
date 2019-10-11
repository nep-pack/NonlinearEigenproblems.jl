using NonlinearEigenproblems, LinearAlgebra

# test on the following DEP
# M(λ)=-λI+A0+exp(-λ)vv'

n=2;
A0=rand(n,n)
v=rand(n,n)


Av=[A0]
Bv=[-one(A0)-v*v']
BvLR=[v/2, v/3, v/4, v/5]
AvLR=[zero(v),zero(v),zero(v),zero(v)]
Z=v;

d=4;
M=diagm( 0 =>  ones(d) )[2:end,:]
N=diagm( -1 =>  1 ./ (1:d-1) )[2:end,:]

c=CORKPencilLR(M,N,Av,AvLR,Bv,BvLR,Z)

# hardcoded contruction of the linearization
dd=length(Av)
M11=M[1:dd,1:dd]; 			N11=N[1:dd,1:dd];
M21=M[dd+1:end,1:dd];		N21=N[dd+1:end,1:dd];
M22=M[dd+1:end,dd+1:end];	N22=N[dd+1:end,dd+1:end];
O=kron(M[1:dd,dd+1:end],zeros(n,d-dd));
II1=one(A0);	II2=[1]; # II2=one(Z)
AA=	[hcat(Av...) hcat(AvLR...);
	[kron(M11,II1) 			zeros(size(kron(M11,II1),1),size(kron(M22,II2),2)); 
	 kron(M21,copy(Z')) 	kron(M22,II2)]]
#BB=[hcat(Bv...) hcat(BvLR...);]

#AA=[hcat(Av...) hcat(AvLR...); kron(M,II)], [hcat(Bv...); kron(N,II)]
