####################################
# Rational Krylov utility functions
####################################

export lejabagby, scgendivdiffs, ratnewtoncoeffs, ratnewtoncoeffsm

"""
Generate Leja-Bagby points (a,b) on (A,B), with scaling factors β such that
the uniform norm on the control set C is 1. Greedy search for a minimum is
performed on B. If keepA is true then the points in the output a will be exactly
those of A, otherwise the points in a are also chosen via greedy search on A. If
forceInf is a positive integer, the first forceInf poles in b will be infinity.
"""
function lejabagby(A::AbstractVector{CT}, B::AbstractVector{T}, C::AbstractVector{CT}, m::Int, keepA::Bool=false, forceInf::Int=0) where {T<:Real, CT<:Complex{T}}
    if minimum(abs.(B)) < 1e-9
        @warn "There is at least one pole candidate in B being nearby zero. Consider shifting your problem for stability."
    end

    a = [A[1]]
    b = [forceInf > 0 ? T(Inf) : B[1]]
    β = [T(1.)]

    o = one(eltype(A))
    sA = fill(o, size(A))
    sB = fill(o, size(B))
    sC = fill(o, size(C))

    @inbounds for j = 1:m-1
        binv = 1 / b[j]
        βinv = 1 / β[j]
        for k = 1:length(A); sA[k] = sA[k] * βinv * (A[k] - a[j]) / (1 - A[k] * binv); end
        for k = 1:length(B); sB[k] = sB[k] * βinv * (B[k] - a[j]) / (1 - B[k] * binv); end
        for k = 1:length(C); sC[k] = sC[k] * βinv * (C[k] - a[j]) / (1 - C[k] * binv); end

        push!(a, A[keepA ? j+1 : argmax([isnan(x) ? -Inf : abs.(x) for x in sA])])
        push!(b, forceInf > j ? Inf : B[argmin([isnan(x) ? Inf : abs.(x) for x in sB])])
        push!(β, maximum(abs.(sC)))

        # treat single point case
        if β[j+1] < eps()
            β[j+1] = 1
        end
    end

    return a, b, β
end

"""
Compute scalar generalized divided differences.

# Arguments
- `σ`: Discretization of target set.
- `ξ`: Discretization of singularity set.
- `β`: Scaling factors.
"""
function scgendivdiffs(σ::AbstractVector{CT}, ξ, β, maxdgr, isfunm, pff) where CT<:Complex{<:Real}
    sgdd = zeros(CT, length(pff), maxdgr+2)
    for ii = 1:length(pff)
        if isfunm
            sgdd[ii,:] = ratnewtoncoeffsm(pff[ii], σ, ξ, β)
        else
            sgdd[ii,:] = map(m -> m[1], ratnewtoncoeffs(pff[ii], σ, ξ, β))
        end
    end
    return sgdd
end

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
