using NonlinearEigenproblems, Random, LinearAlgebra

# representation of the structured linearizations used in the CORK framework
struct CORK_pencil{T<:AbstractMatrix}
    M::T
    N::T
    Av::Vector{T}   # Array of Array of matrices
    Bv::Vector{T}   # Array of Array of matrices
end

n=2
M=rand(n,n)
N=rand(n,n)
kk=3
Av=[rand(3,3) rand(3,3)]
Bv=[rand(3,3) rand(3,3)]
