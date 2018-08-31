# Gun: variant R2 (fully rational case; leja + repeated nodes + freeze)

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

@testset "NLEIGS: Gun variant R2" begin
    verbose = 1

    nep, Σ, Ξ, v, nodes, funres = gun_init()

    # solve nlep
    @time lambda, X, res, solution_info = nleigs(nep, Σ, Ξ=Ξ, displaylevel=verbose > 0 ? 1 : 0, minit=60, maxit=100, v=v, nodes=nodes, errmeasure=funres, return_details=verbose > 1)

    nleigs_verify_lambdas(21, nep, X, lambda)

    if verbose > 1
        include("nleigs_residual_plot.jl")
        nleigs_residual_plot("Gun: variant R2", solution_info, Σ; ylims=[1e-17, 1e-1])
    end
end
