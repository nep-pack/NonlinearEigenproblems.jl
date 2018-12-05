using SparseArrays, Random, LinearAlgebra
n=20

K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A = sparse(K, J, rand(3*n-2));

Af=factorize(A)

O=zero(A)
B=[ O A; A' O ]

Bf=factorize(B)
