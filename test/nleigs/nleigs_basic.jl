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
B = Array{Array{Float64,2}}([[1 3; 5 6], [3 4; 6 6], [1 0; 0 1]])
NLEP = Dict("B" => B)

funres = (λ, X) -> map(i -> norm(B[1]*X[:,i] + λ[i]*(B[2]*X[:,i]) + λ[i]^2*(B[3]*X[:,i])), 1:length(λ))

Sigma = [-10.0-2im, 10-2im, 10+2im, -10+2im]

@testset "NLEIGS: Polynomial only" begin
    options = Dict("maxit" => 10, "v0" => ones(n), "funres" => funres)
    @time X, lambda = nleigs(NLEP, Sigma, options=options)
    nleigs_verify_lambdas(4, NLEP, X, lambda)
end

@testset "NLEIGS: Non-convergent linearization" begin
    options = Dict("maxit" => 10, "v0" => ones(n), "maxdgr" => 5, "funres" => funres)
    @time X, lambda = nleigs(NLEP, Sigma, options=options)
    nleigs_verify_lambdas(4, NLEP, X, lambda)
end

@testset "NLEIGS: Non-convergent linearization (static)" begin
    options = Dict("maxit" => 10, "v0" => ones(n), "maxdgr" => 5, "funres" => funres, "static" => true)
    @time X, lambda = nleigs(NLEP, Sigma, options=options)
    nleigs_verify_lambdas(4, NLEP, X, lambda)
end

@testset "NLEIGS: return_info" begin
    options = Dict("maxit" => 10, "v0" => ones(n), "funres" => funres)
    @time X, lambda, res, solution_info = nleigs(NLEP, Sigma, options=options, return_info=true)
    nleigs_verify_lambdas(4, NLEP, X, lambda)

    info_λ = solution_info["Lam"][:,end]
    in_sigma = map(p -> inpolygon(real(p), imag(p), real(Sigma), imag(Sigma)), info_λ)
    info_λ = info_λ[in_sigma]

    # test that eigenvalues in the info are the same as those returned by nleigs
    @test length(length(info_λ)) == 4
    @test length(union(lambda, info_λ)) == 4

    # test that the residuals are near 0
    info_res = solution_info["Res"][in_sigma,end]
    @test all(r -> r < 1e-12, info_res)
end
