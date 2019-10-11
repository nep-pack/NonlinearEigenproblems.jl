using NonlinearEigenproblems

# test on the following DEP
# M(λ)=-λI+A0+exp(-λ)vv'

A0=rand(2,2)
v=rand(2,1)


Av=[A0]
Bv=[-one(A0)-v*v']
BvRL=[v/2, v/3, v/4]
AvRL=[zero(v),zero(v),zero(v)]
Z=copy(v');
M=rand(2,2)
N=rand(2,2)

c=CORKPencilLR(M,N,Av,AvLR,Bv,BvLR,Z)
