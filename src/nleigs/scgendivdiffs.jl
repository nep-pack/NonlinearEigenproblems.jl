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
