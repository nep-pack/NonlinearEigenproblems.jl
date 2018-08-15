# Particle: variant S

# Intended to be run from nep-pack/ directory or nep-pack/test directory
if !isdefined(:global_modules_loaded)
    workspace()

    push!(LOAD_PATH, string(@__DIR__, "/../../src"))

    using NEPCore
    using NEPTypes
    using Gallery
    using Base.Test
end

include("nleigs_test_utils.jl")
include("particle_test_utils.jl")
include("../../src/nleigs/method_nleigs.jl")

verbose = 1

nep, Sigma, Xi, v0, nodes, xmin, xmax = particle_init(2)

# solve nlep
@time X, lambda, res, solution_info = nleigs(nep, Sigma, Xi=Xi, verbose=verbose > 0 ? 1 : 0, maxdgr=50, minit=120, maxit=200, v0=v0, nodes=nodes, static=true, return_details=verbose > 1)

@testset "NLEIGS: Particle variant S" begin
    nleigs_verify_lambdas(2, nep, X, lambda)
end

if verbose > 1
    include("nleigs_residual_plot.jl")
    approxSigma = [xmin-im*1e-10, xmin+im*1e-10, xmax+im*1e-10, xmax-im*1e-10]
    nleigs_residual_plot("Particle: variant S", solution_info, approxSigma)
end
