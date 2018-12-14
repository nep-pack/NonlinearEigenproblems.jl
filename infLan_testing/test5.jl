using SparseArrays, Random, LinearAlgebra, BenchmarkTools
n=20

K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A = sparse(K, J, rand(3*n-2));

Af=factorize(A)

O=zero(A)
M=[ O A; A' O ]

try
    MM=factorize(M)
catch e
    println("t ",typeof(e))
    if PosDefException
        MM=lu(M)
    end
end
b=rand(size(M,1))
@btime begin x=M\b; end
@btime begin xx=MM\b; end
#norm(M*xx-b)
