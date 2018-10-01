using NonlinearEigenproblems
using NonlinearEigenproblems.RKHelper

"""
Construct a pencil for the linearization of a nonlinear eigenvalue problem.

# Arguments
- `nep`: An instance of a nonlinear eigenvalue problem.
- `Σ`: A vector containing the points of a polygonal target set in the complex plane.
- `Ξ`: A vector containing a discretization of the singularity set.
- `maxdgr`: Max degree of approximation.
- `tollin`: Tolerance for convergence of linearization.
"""
function linearize(
        ::Type{T},
        nep::NEP,
        Σ::AbstractVector{CT},
        Ξ::Vector{T},
        maxdgr::Int,
        tollin::T) where {T<:Real, CT<:Complex{T}}

    P = get_rk_nep(T, nep)

    # Discretization of Σ --> Gamma & Leja-Bagby points
    gamma,_ = discretizepolygon(Σ)
    σ,ξ,β = lejabagby(gamma, Ξ, gamma, maxdgr+2, false, P.p)

    # Rational Newton coefficients -- compute scalar generalized divided differences
    drange = 1:maxdgr+2
    sgdd = scgendivdiffs(σ[drange], ξ[drange], β[drange], maxdgr, true, get_fv(nep))

    n = size(nep, 1)
    BC = get_Av(nep)
    D = Vector()
    Dt = Vector()

    CONVERGENCE_WINDOW = 5
    norms = zeros(T, CONVERGENCE_WINDOW, 2)

    for nb = 0:maxdgr
        # full rank
        Dk = spzeros(CT, n, n)
        for ii = 1:(P.p+1+P.q)
            Dk += sgdd[ii,nb+1] * BC[ii]
        end
        push!(D, Dk)

        # low rank
        if nb > P.p
            for ii = 1:P.q
                d = sgdd[P.p+1+ii,nb+1] * P.L[ii]
                Dtk = ii == 1 ? d : hcat(Dtk, d)
            end
        else
            Dtk = Matrix{CT}(undef,0,0)
        end
        push!(Dt, Dtk)

        # break when the CONVERGENCE_WINDOW most recent norms are below the tolerance
        idx = 1 + mod(nb, size(norms, 1))
        norms[idx, :] = [norm(Dk), norm(Dtk)]
        all(mean(norms, dims=1) .< tollin) && break
    end

    d = length(D) - 1
    return ConvergedLinearization(n, P.p, P.UU, d, σ[1:d+1], ξ[1:d+1], β[1:d+1], D, Dt)
end
