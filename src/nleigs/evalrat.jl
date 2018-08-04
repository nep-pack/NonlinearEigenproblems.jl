#EVALRAT Evaluate nodal rational function at the points z
function evalrat(sigma, xi, beta, z::AbstractVector{<:Number})
    r = ones(size(z)) / beta[1]
    for j = 1:length(sigma)
        r .*= (z - sigma[j]) ./ (1 - z/xi[j]) / beta[j+1]
    end
    return r
end
