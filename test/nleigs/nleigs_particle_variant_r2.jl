# Particle: variant R2 (fully rational case; leja + repeated nodes + freeze)

push!(LOAD_PATH, normpath(@__DIR__, "..")); using TestUtils
using NonlinearEigenproblems
using Test

include("nleigs_test_utils.jl")
include("particle_test_utils.jl")

@bench @testset "NLEIGS: Particle variant R2" begin
    verbose = displaylevel

    nep, Σ, Ξ, v, nodes, xmin, xmax = particle_init(2)

    # solve nlep
    lambda, X, res, solution_info = nleigs(nep, Σ, Ξ=Ξ, displaylevel=verbose > 0 ? 1 : 0, maxdgr=50, minit=30, maxit=100, v=v, nodes=nodes, return_details=verbose > 1)

    nleigs_verify_lambdas(2, nep, X, lambda)

    if verbose > 1
        include("nleigs_residual_plot.jl")
        approx_Σ = [xmin-im*1e-10, xmin+im*1e-10, xmax+im*1e-10, xmax-im*1e-10]
        nleigs_residual_plot("Particle: variant R2", solution_info, approx_Σ)
    end
end
