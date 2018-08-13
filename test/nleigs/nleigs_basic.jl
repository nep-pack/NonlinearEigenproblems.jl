# Solves a few basic eigenvalue problems to test various aspects of NLEIGS

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
include("../../src/nleigs/method_nleigs.jl")

n = 2
B = Vector{Matrix{Float64}}([[1 3; 5 6], [3 4; 6 6], [1 0; 0 1]])
nep = SPMFLowRankNEP(B, Vector{SPMFLowRankMatrix{Matrix{Float64}}}(0))

funres = (λ, X) -> map(i -> norm(B[1]*X[:,i] + λ[i]*(B[2]*X[:,i]) + λ[i]^2*(B[3]*X[:,i])), 1:length(λ))

Sigma = [-10.0-2im, 10-2im, 10+2im, -10+2im]

@testset "NLEIGS: Polynomial only" begin
    options = Dict("maxit" => 10, "v0" => ones(n), "funres" => funres, "blksize" => 5)
    @time X, lambda = nleigs(nep, Sigma, options=options)
    nleigs_verify_lambdas(4, nep, X, lambda)
end

@testset "NLEIGS: Non-convergent linearization" begin
    options = Dict("maxit" => 10, "v0" => ones(n), "maxdgr" => 5, "funres" => funres, "blksize" => 5)
    @time X, lambda = nleigs(nep, Sigma, options=options)
    nleigs_verify_lambdas(4, nep, X, lambda)
end

@testset "NLEIGS: Non-convergent linearization (static)" begin
    options = Dict("maxit" => 10, "v0" => ones(n), "maxdgr" => 5, "funres" => funres, "static" => true, "blksize" => 5)
    @time X, lambda = nleigs(nep, Sigma, options=options)
    nleigs_verify_lambdas(4, nep, X, lambda)
end

@testset "NLEIGS: return_info" begin
    options = Dict("maxit" => 10, "v0" => ones(n), "funres" => funres, "blksize" => 5)
    @time X, lambda, res, solution_info = nleigs(nep, Sigma, options=options, return_info=true)
    nleigs_verify_lambdas(4, nep, X, lambda)

    info_λ = solution_info["Lam"][:,end]
    local in_sigma = map(p -> inpolygon(real(p), imag(p), real(Sigma), imag(Sigma)), info_λ)
    info_λ = info_λ[in_sigma]

    # test that eigenvalues in the info are the same as those returned by nleigs
    @test length(info_λ) == 4
    @test length(union(lambda, info_λ)) == 4

    # test that the residuals are near 0
    info_res = solution_info["Res"][in_sigma,end]
    @test all(r -> r < 1e-12, info_res)
end

struct TestNEP <: NEP
    n::Int  # problem size; this is the only required field in the custom NEP type
end

# compute_Mder (for the 0:th derivative) has to be implemented to solve a custom NEP type with NLEIGS
import NEPCore.compute_Mder
compute_Mder(_::TestNEP, λ::Number) = compute_Mder(nep, λ)

# compute_Mlincomb is needed for the test verification only
import NEPCore.compute_Mlincomb
compute_Mlincomb(_::TestNEP, λ::Number, v) = compute_Mlincomb(nep, λ, v)

@testset "NLEIGS: Custom function" begin
    testnep = TestNEP(nep.n)

    options = Dict("maxit" => 10, "v0" => ones(n), "funres" => funres, "blksize" => 5)
    @time X, lambda = nleigs(testnep, Sigma, options=options)
    nleigs_verify_lambdas(4, nep, X, lambda)
end
