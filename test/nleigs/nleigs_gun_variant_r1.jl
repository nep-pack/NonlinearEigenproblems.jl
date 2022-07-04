# Gun: variant R1 (fully rational case; only repeated nodes)

#using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test

include(joinpath("..", "rk_helper", "gun_test_utils.jl"))

@bench @testset "NLEIGS: Gun variant R1" begin
    verbose = displaylevel

    nep, Σ, Ξ, v, nodes, funres = gun_init()

    # solve nlep
    lambda, X, res, solution_info = nleigs(nep, Σ, Ξ=Ξ, logger=verbose > 0 ? 1 : 0, maxit=100, v=v, leja=0, nodes=nodes, reusefact=2, errmeasure=funres, return_details=verbose > 1)

    verify_lambdas(21, nep, lambda, X)

    if verbose > 1
        include("nleigs_residual_plot.jl")
        nleigs_residual_plot("Gun: variant R1", solution_info, Σ; ylims=[1e-17, 1e-1])
    end
end
