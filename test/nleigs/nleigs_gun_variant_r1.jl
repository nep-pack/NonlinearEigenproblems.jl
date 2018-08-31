# Gun: variant R1 (fully rational case; only repeated nodes)

# Intended to be run from nep-pack/ directory or nep-pack/test directory
push!(LOAD_PATH, string(@__DIR__, "/../../src"))

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver
using Gallery
using IterativeSolvers
using Test

include("nleigs_test_utils.jl")
include("gun_test_utils.jl")

@testset "NLEIGS: Gun variant R1" begin
    verbose = 1

    nep, Σ, Ξ, v, nodes, funres = gun_init()

    # solve nlep
    @time lambda, X, res, solution_info = nleigs(nep, Σ, Ξ=Ξ, displaylevel=verbose > 0 ? 1 : 0, maxit=100, v=v, leja=0, nodes=nodes, reuselu=2, errmeasure=funres, return_details=verbose > 1)

    nleigs_verify_lambdas(21, nep, X, lambda)

    if verbose > 1
        include("nleigs_residual_plot.jl")
        nleigs_residual_plot("Gun: variant R1", solution_info, Σ; ylims=[1e-17, 1e-1])
    end
end
