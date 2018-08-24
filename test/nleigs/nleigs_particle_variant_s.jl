# Particle: variant S

# Intended to be run from nep-pack/ directory or nep-pack/test directory
if !isdefined(:global_modules_loaded)
    workspace()

    push!(LOAD_PATH, string(@__DIR__, "/../../src"))
    push!(LOAD_PATH, string(@__DIR__, "/../../src/nleigs"))

    using NEPCore
    using NEPTypes
    using NleigsTypes
    using Gallery
    using IterativeSolvers
    using Base.Test
end

include("nleigs_test_utils.jl")
include("particle_test_utils.jl")
include("../../src/nleigs/method_nleigs.jl")

verbose = 1

nep, Σ, Ξ, v, nodes, xmin, xmax = particle_init(2)

# solve nlep
@time X, lambda, res, solution_info = nleigs(nep, Σ, Ξ=Ξ, displaylevel=verbose > 0 ? 1 : 0, maxdgr=50, minit=120, maxit=200, v=v, nodes=nodes, static=true, return_details=verbose > 1)

@testset "NLEIGS: Particle variant S" begin
    nleigs_verify_lambdas(2, nep, X, lambda)
end

if verbose > 1
    include("nleigs_residual_plot.jl")
    approx_Σ = [xmin-im*1e-10, xmin+im*1e-10, xmax+im*1e-10, xmax-im*1e-10]
    nleigs_residual_plot("Particle: variant S", solution_info, approx_Σ)
end
