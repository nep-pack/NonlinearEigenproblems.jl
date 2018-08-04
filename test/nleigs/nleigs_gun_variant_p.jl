# Gun: variant P (polynomial case; only repeated nodes)

using Base.Test

include("nleigs_test_utils.jl")
include("gun_test_utils.jl")
include("../../src/nleigs/method_nleigs.jl")

verbose = false

NLEP, Sigma, Xi, v0, nodes, funres = gun_init()

Xi = Float64[] # no pole candidates

options = Dict(
    "disp" => verbose ? 1 : 0,
    "maxit" => 100,
    "v0" => v0,
    "funres" => funres,
    "leja" => 0,
    "nodes" => nodes,
    "reuselu" => 2)

# solve nlep
@time X, lambda, res, solution_info = nleigs(NLEP, Sigma, Xi=Xi, options=options, return_info=verbose)

@testset "NLEIGS: Gun variant P" begin
    nleigs_verify_lambdas(17, NLEP, X, lambda)
end

if verbose
    include("nleigs_residual_plot.jl")
    nleigs_residual_plot("Gun: variant P", solution_info, Sigma; ylims=[1e-17, 1e-1])
end
