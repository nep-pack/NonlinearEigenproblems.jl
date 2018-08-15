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
nep = PNEP(B, Vector{MatrixAndFunction{Matrix{Float64}}}(0))

Sigma = [-10.0-2im, 10-2im, 10+2im, -10+2im]

@testset "NLEIGS: Polynomial only" begin
    options = Dict("maxit" => 10, "v0" => ones(n), "blksize" => 5)
    @time X, lambda = nleigs(nep, Sigma, options=options)
    nleigs_verify_lambdas(4, nep, X, lambda)
end

@testset "NLEIGS: Non-convergent linearization" begin
    options = Dict("maxit" => 10, "v0" => ones(n), "maxdgr" => 5, "blksize" => 5)
    @test_warn "Linearization not converged" begin
        @time X, lambda = nleigs(nep, Sigma, options=options)
        nleigs_verify_lambdas(4, nep, X, lambda)
    end
end

@testset "NLEIGS: Non-convergent linearization (static)" begin
    options = Dict("maxit" => 10, "v0" => ones(n), "maxdgr" => 5, "static" => true, "blksize" => 5)
    @test_warn "Linearization not converged" begin
        @time X, lambda = nleigs(nep, Sigma, options=options)
        nleigs_verify_lambdas(4, nep, X, lambda)
    end
end

@testset "NLEIGS: return_details" begin
    options = Dict("maxit" => 10, "v0" => ones(n), "blksize" => 5)
    @time X, lambda, res, details = nleigs(nep, Sigma, options=options, return_details=true)
    nleigs_verify_lambdas(4, nep, X, lambda)

    info_λ = details.Lam[:,end]
    local in_sigma = map(p -> inpolygon(real(p), imag(p), real(Sigma), imag(Sigma)), info_λ)
    info_λ = info_λ[in_sigma]

    # test that eigenvalues in the info are the same as those returned by nleigs
    @test length(info_λ) == 4
    @test length(union(lambda, info_λ)) == 4

    # test that the residuals are near 0
    info_res = details.Res[in_sigma,end]
    @test all(r -> r < 1e-12, info_res)
end
