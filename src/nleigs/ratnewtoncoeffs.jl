include("evalrat.jl")

"""
Compute rational divided differences for the function fun (can be matrix
valued), using differencing. The sigmas need to be distinct. For scalar
functions or non-distinct sigmas it may be better to use `ratnewtoncoeffsm`.
"""
function ratnewtoncoeffs(fun, sigma::AbstractVector{CT}, xi::AbstractVector{T}, beta::AbstractVector{T}) where {T<:Real, CT<:Complex{T}}
    m = length(sigma)
    D = Vector{Matrix{CT}}(m)

    # compute divided differences D0,D1,...,Dm
    as_matrix(x::Number) = (M = Matrix{eltype(x)}(1,1); M[1] = x; M)
    D[1] = fun(as_matrix(sigma[1])) * beta[1]
    n = size(D[1], 1)
    for j = 2:m
        # evaluate current linearizaion at sigma[j]
        Qj = T(0)
        for k = 1:j-1
            Qj += D[k] * evalrat(sigma[1:k-1], xi[1:k-1], beta[1:k], [sigma[j]])[1]
        end

        # get divided difference from recursion (could be done via Horner)
        D[j] = (fun(as_matrix(sigma[j])) - Qj) / evalrat(sigma[1:j-1], xi[1:j-1], beta[1:j], [sigma[j]])[1]
    end

    return D
end
