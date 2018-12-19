using LinearAlgebra, ToeplitzMatrices, BenchmarkTools

# construct a random Hankel matrix
m=100
n=100
p=20
T=Float64
A=Vector{Hankel{Float64}}(undef,p)
for j=1:p
    v=rand(m)
    w=rand(m)
    v[end]=w[1]
    A[j]=Hankel(v,w)
end
