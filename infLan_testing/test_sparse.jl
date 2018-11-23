using Random, SparseArrays, BenchmarkTools

n=100000
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A = sparse(K, J, rand(3*n-2)); A = A+A';

v1=rand(n,1)
v2=rand(n,1)
v3=rand(n,1)

V=[v1 v2 v3]

@btime begin
    z1=A*v1
    z2=A*v2
    z3=A*v3
end

@btime begin
    Z=A*V
end
1
