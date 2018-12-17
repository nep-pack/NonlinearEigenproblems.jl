using LinearAlgebra, ToeplitzMatrices

# construct a random Hankel matrix
v=rand(n)
w=rand(n)
v[end]=w[1]
A=Hankel(v,w)
X=rand(size(A,2))
Y1=A*X

AA=Matrix(A)
Y2=AA*X

norm(Y1-Y2)
