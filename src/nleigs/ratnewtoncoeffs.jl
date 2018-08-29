include("evalrat.jl")

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
        Qj = T(0)
        for k = 1:j-1
            Qj += D[k][1] * evalrat(σ[1:k-1], ξ[1:k-1], β[1:k], [σ[j]])[1]
        end

        # get divided difference from recursion (could be done via Horner)
        D[j] = (fun(as_matrix(σ[j])) .- Qj) / evalrat(σ[1:j-1], ξ[1:j-1], β[1:j], [σ[j]])[1]
    end

    return D
end
