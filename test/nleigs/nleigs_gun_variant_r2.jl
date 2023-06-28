# Gun: variant R2 (fully rational case; leja + repeated nodes + freeze)

#using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test

include(joinpath("..", "rk_helper", "gun_test_utils.jl"))

@bench @testset "NLEIGS: Gun variant R2" begin
    verbose = displaylevel

    nep, Σ, Ξ, v, nodes, funres = gun_init()

    # solve nlep
    lambda, X, res, solution_info = nleigs(nep, Σ, Ξ=Ξ, logger=verbose > 0 ? 1 : 0, minit=60, maxit=100, v=v, nodes=nodes, errmeasure=funres, return_details=verbose > 1)

    verify_lambdas(21, nep, lambda, X)

    if verbose > 1
        include("nleigs_residual_plot.jl")
        nleigs_residual_plot("Gun: variant R2", solution_info, Σ; ylims=[1e-17, 1e-1])
    end
end
