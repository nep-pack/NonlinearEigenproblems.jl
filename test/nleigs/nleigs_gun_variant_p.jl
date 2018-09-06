# Gun: variant P (polynomial case; only repeated nodes)

using NonlinearEigenproblems.NEPSolver
using Test

include("nleigs_test_utils.jl")
include("gun_test_utils.jl")

@testset "NLEIGS: Gun variant P" begin
    verbose = 1

    nep, Σ, _, v, nodes, funres = gun_init()

    # solve nlep
    @time lambda, X, res, solution_info = nleigs(nep, Σ, displaylevel=verbose > 0 ? 1 : 0, maxit=100, v=v, leja=0, nodes=nodes, reuselu=2, errmeasure=funres, return_details=verbose > 1)

    nleigs_verify_lambdas(17, nep, X, lambda)

    if verbose > 1
        include("nleigs_residual_plot.jl")
        nleigs_residual_plot("Gun: variant P", solution_info, Σ; ylims=[1e-17, 1e-1])
    end
end
