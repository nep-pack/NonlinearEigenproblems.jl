#EVALRAT Evaluate nodal rational function at the points z
function evalrat(sigma::AbstractVector{Complex{T}}, xi::AbstractVector{T}, beta::AbstractVector{T}, z::AbstractVector{Complex{T}}) where T<:Real
    r = complex(ones(size(z)) / beta[1])
    for j = 1:length(sigma)
        r .*= (z - sigma[j]) ./ (1 - z/xi[j]) / beta[j+1]
    end
    return r
end
