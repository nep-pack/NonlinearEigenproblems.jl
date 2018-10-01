"""
Compute rational divided differences for the function fun (can be matrix
valued), using differencing. The σ need to be distinct. For scalar
functions or non-distinct σ it may be better to use `ratnewtoncoeffsm`.
"""
function ratnewtoncoeffs(fun, σ::AbstractVector{CT}, ξ::AbstractVector{T}, β::AbstractVector{T}) where {T<:Real, CT<:Complex{T}}
    m = length(σ)
    D = Vector{Matrix{CT}}(undef, m)

    # compute divided differences D0,D1,...,Dm
    as_matrix(x::Number) = (M = Matrix{eltype(x)}(undef,1,1); M[1] = x; M)
    D[1] = fun(as_matrix(σ[1])) * β[1]

    for j = 2:m
        # evaluate current linearizaion at σ[j]
        Qj = zeros(T, size(D[1]))
        for k = 1:j-1
            Qj += D[k] * evalrat(σ[1:k-1], ξ[1:k-1], β[1:k], [σ[j]])[1]
        end

        # get divided difference from recursion (could be done via Horner)
        D[j] = (fun(as_matrix(σ[j])) .- Qj) / evalrat(σ[1:j-1], ξ[1:j-1], β[1:j], [σ[j]])[1]
    end

    return D
end

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

"Evaluate nodal rational function at the points z."
function evalrat(σ::AbstractVector{CT}, ξ::AbstractVector{T}, β::AbstractVector{T}, z::AbstractVector{CT}) where {T<:Real, CT<:Complex{T}}
    r = ones(CT, size(z)) / β[1]
    for j = 1:length(σ)
        r .*= (z .- σ[j]) ./ (1 .- z/ξ[j]) / β[j+1]
    end
    return r
end
