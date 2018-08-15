# Compute rational divided differences for the SCALAR function fun, using
# matrix functions. fm has to be a handle to a matrix function fm(A).
function ratnewtoncoeffsm(fm, sigma::AbstractVector{T}, xi::AbstractVector{RT}, beta::AbstractVector{RT}) where {T<:Number, RT<:Real}
    m = length(sigma) - 1

    sigma = sigma[:]
    xi = xi[:]
    beta = beta[:]

    # build Hessenberg matrices
    K = Bidiagonal(ones(m+1), beta[2:m+1]./xi[1:m], 'L')
    H = Bidiagonal(sigma[1:m+1], beta[2:m+1], 'L')

    # column balancing
    P = Diagonal(1./maximum(abs.(K), 1)[:])
    K *= P
    H *= P

    D = fm(H/K) * eye(m+1, 1) * beta[1]
    D = D.'

    return D
end
