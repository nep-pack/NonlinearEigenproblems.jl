using NonlinearEigenproblems

include("../nleigs/method_nleigs.jl")

"""
Constructs a pencil for the linearization of a nonlinear eigenvalue problem.

# Arguments
- `nep`: An instance of a nonlinear eigenvalue problem.
- `Σ`: A vector containing the points of a polygonal target set in the complex plane.
- `Ξ`: A vector containing a discretization of the singularity set.
- `maxdgr`: Max degree of approximation.
- `tollin`: Tolerance for convergence of linearization.
"""
function linearize(
        ::Type{T},
        nep::NleigsNEP,
        Σ::AbstractVector{CT},
        Ξ::Vector{T},
        maxdgr::Int,
        tollin::T) where {T<:Real, CT<:Complex{T}}

    # Discretization of Σ --> Gamma & Leja-Bagby points
    gamma,_ = discretizepolygon(Σ)
    σ,ξ,β = lejabagby(gamma, Ξ, gamma, maxdgr+2, false, nep.p)

    # Rational Newton coefficients -- compute scalar generalized divided differences
    drange = 1:maxdgr+2
    sgdd = scgendivdiffs(σ[drange], ξ[drange], β[drange], maxdgr, true, get_fv(nep.nep))

    BC = get_Av(P.nep)

    n = size(nep.nep, 1)
    nrmD = Vector()
    nrmDt = Vector()
    D = Vector()
    Dt = Vector()

    for nb = 0:maxdgr
        Dk = spzeros(CT, n, n)
        for ii = 1:(nep.p+1+nep.q)
            Dk += sgdd[ii,nb+1] * BC[ii]
        end
        push!(D, Dk)
        push!(nrmD, norm(Dk))

        if nb > nep.p
            for ii = 1:nep.q
                d = sgdd[nep.p+1+ii,nb+1] * nep.L[ii]
                Dtk = ii == 1 ? d : hcat(Dtk, d)
            end
        else
            Dtk = Matrix{CT}(undef,0,0)
        end
        push!(Dt, Dtk)
        push!(nrmDt, norm(Dtk))

        if sum(nrmD[max(1,end-4):end]) < 5 * tollin && sum(nrmDt[max(1,end-4):end]) < 5 * tollin
            break
        end
    end

    d = length(D) - 1
    return Linearization(n, nep.p, nep.UU, d, σ[1:d+1], ξ[1:d+1], β[1:d+1], D, Dt)
end

struct Linearization
    n::Int # Size of problem
    p::Int # Order of polynomial part
    U::AbstractMatrix # U factors of the low rank nonlinear part
    d::Int # Degree of approximation
    σ::Vector # Interpolation nodes
    ξ::Vector # Poles
    β::Vector # Scaling factors
    D::Vector #{<:AbstractMatrix} # Full rank generalized divided differences
    Dt::Vector # Low rank generalized divided differences
end
