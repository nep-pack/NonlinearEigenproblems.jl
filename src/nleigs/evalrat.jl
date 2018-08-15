#EVALRAT Evaluate nodal rational function at the points z
function evalrat(sigma::AbstractVector{T}, xi::AbstractVector{RT}, beta::AbstractVector{RT}, z::AbstractVector{T}) where {T<:Number, RT<:Real}
    r = ones(T, size(z)) / beta[1]
    for j = 1:length(sigma)
        r .*= (z - sigma[j]) ./ (1 - z/xi[j]) / beta[j+1]
    end
    return r
end
