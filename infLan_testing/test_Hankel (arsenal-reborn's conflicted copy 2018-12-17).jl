using LinearAlgebra, ToeplitzMatrices, BenchmarkTools

# construct a random Hankel matrix
m=100
n=100
v=rand(m)
w=rand(m)
v[end]=w[1]
A=Hankel(v,w)
X=rand(n,size(A,2))
@btime begin Y1=X*A end

AA=Matrix(A)
@btime begin Y2=X*AA end

norm(Y1-Y2)
