"Evaluate nodal rational function at the points z."
function evalrat(σ::AbstractVector{CT}, ξ::AbstractVector{T}, β::AbstractVector{T}, z::AbstractVector{CT}) where {T<:Real, CT<:Complex{T}}
    r = ones(CT, size(z)) / β[1]
    for j = 1:length(σ)
        r .*= (z .- σ[j]) ./ (1 .- z/ξ[j]) / β[j+1]
    end
    return r
end
