"Evaluate nodal rational function at the points z."
function evalrat(sigma::AbstractVector{CT}, xi::AbstractVector{T}, beta::AbstractVector{T}, z::AbstractVector{CT}) where {T<:Real, CT<:Complex{T}}
    r = ones(CT, size(z)) / beta[1]
    for j = 1:length(sigma)
        r .*= (z - sigma[j]) ./ (1 - z/xi[j]) / beta[j+1]
    end
    return r
end
