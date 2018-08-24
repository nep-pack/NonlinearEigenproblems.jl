"Evaluate nodal rational function at the points z."
function evalrat(σ::AbstractVector{CT}, xi::AbstractVector{T}, beta::AbstractVector{T}, z::AbstractVector{CT}) where {T<:Real, CT<:Complex{T}}
    r = ones(CT, size(z)) / beta[1]
    for j = 1:length(σ)
        r .*= (z - σ[j]) ./ (1 - z/xi[j]) / beta[j+1]
    end
    return r
end
