# Particle: variant S

using NonlinearEigenproblems.NEPSolver
using Test

include("nleigs_test_utils.jl")
include("particle_test_utils.jl")

@testset "NLEIGS: Particle variant S" begin
    verbose = 1

    nep, Σ, Ξ, v, nodes, xmin, xmax = particle_init(2)

    # solve nlep
    @time lambda, X, res, solution_info = nleigs(nep, Σ, Ξ=Ξ, displaylevel=verbose > 0 ? 1 : 0, maxdgr=50, minit=120, maxit=200, v=v, nodes=nodes, static=true, return_details=verbose > 1)

    nleigs_verify_lambdas(2, nep, X, lambda)

    if verbose > 1
        include("nleigs_residual_plot.jl")
        approx_Σ = [xmin-im*1e-10, xmin+im*1e-10, xmax+im*1e-10, xmax-im*1e-10]
        nleigs_residual_plot("Particle: variant S", solution_info, approx_Σ)
    end
end
