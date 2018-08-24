"""
Compute rational divided differences for the scalar function fun, using matrix
functions. `fm` has to be a function object representing a matrix function.
"""
function ratnewtoncoeffsm(fm, σ::AbstractVector{CT}, xi::AbstractVector{T}, β::AbstractVector{T}) where {T<:Real, CT<:Complex{T}}
    m = length(σ) - 1

    σ = σ[:]
    xi = xi[:]
    β = β[:]

    # build Hessenberg matrices
    K = Bidiagonal(ones(m+1), β[2:m+1]./xi[1:m], 'L')
    H = Bidiagonal(σ[1:m+1], β[2:m+1], 'L')

    # column balancing
    P = Diagonal(1./maximum(abs.(K), 1)[:])
    K *= P
    H *= P

    D = fm(H/K) * eye(m+1, 1) * β[1]
    D = D.'

    return D
end
