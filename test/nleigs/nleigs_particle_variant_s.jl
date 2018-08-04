# Particle: variant S

using Base.Test

include("nleigs_test_utils.jl")
include("particle_test_utils.jl")
include("../../src/nleigs/method_nleigs.jl")

verbose = false

NLEP, Sigma, Xi, v0, nodes, funres, xmin, xmax = particle_init(2)

options = Dict(
    "disp" => verbose ? 1 : 0,
    "maxdgr" => 50,
    "minit" => 120,
    "maxit" => 200,
    "v0" => v0,
    "funres" => funres,
    "nodes" => nodes,
    "static" => true)

# solve nlep
@time X, lambda, res, solution_info = nleigs(NLEP, Sigma, Xi=Xi, options=options, return_info=verbose)

@testset "NLEIGS: Particle variant S" begin
    nleigs_verify_lambdas(2, NLEP, X, lambda)
end

if verbose
    include("nleigs_residual_plot.jl")
    approxSigma = [xmin-im*1e-10, xmin+im*1e-10, xmax+im*1e-10, xmax-im*1e-10]
    nleigs_residual_plot("Particle: variant S", solution_info, approxSigma)
end
