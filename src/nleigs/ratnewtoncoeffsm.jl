"""
Compute rational divided differences for the scalar function fun, using matrix
functions. `fm` has to be a function object representing a matrix function.
"""
function ratnewtoncoeffsm(fm, σ::AbstractVector{CT}, ξ::AbstractVector{T}, β::AbstractVector{T}) where {T<:Real, CT<:Complex{T}}
    m = length(σ) - 1

    σ = σ[:]
    ξ = ξ[:]
    β = β[:]

    # build Hessenberg matrices
    K = Bidiagonal(ones(m+1), β[2:m+1]./ξ[1:m], :L)
    H = Bidiagonal(σ[1:m+1], CT.(β[2:m+1]), :L)

    # column balancing
    P = Diagonal(1 ./ maximum(abs.(K), dims = 1)[:])
    K *= P
    H *= P
    H=Matrix(H) # this matrix is sparse, it is needed to convert before matrix-fun evaluation
    D = fm(H/CT.(K))[:,1] * β[1] 
    D = copy(transpose(D))

    return D
end
