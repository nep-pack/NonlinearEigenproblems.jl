# Solves a few basic eigenvalue problems with different NEP types through NLEIGS

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
include("../../src/nleigs/method_nleigs.jl")

n = 2
B = Vector{Matrix{Float64}}([[1 3; 5 6], [3 4; 6 6]])
C = Vector{Matrix{Float64}}([[1 0; 0 1]])
f = [λ -> λ^2]

Sigma = [-10.0-2im, 10-2im, 10+2im, -10+2im]

pep = PEP(B)
spmf = LowRankFactorizedNEP([MatrixAndFunction(C[1], f[1])])
sumnep = SumNEP(pep, spmf)

@testset "NLEIGS: SumNEP" begin
    @time X, lambda = nleigs(sumnep, Sigma, maxit=10, v=ones(n), blksize=5)
    nleigs_verify_lambdas(4, sumnep, X, lambda)
end

@testset "NLEIGS: SPMF_NEP" begin
    spmf_nep = SPMF_NEP([B; C], [λ -> 1; λ -> λ; λ -> λ^2])
    @time X, lambda = nleigs(spmf_nep, Sigma, maxit=10, v=ones(n), blksize=5)
    nleigs_verify_lambdas(4, spmf_nep, X, lambda)
end

struct CustomNLEIGSNEP <: NEP
    n::Int  # problem size; this is the only required field in a custom NEP type when used with NLEIGS
end

# implement a few methods used by the solver
import NEPCore.compute_Mder, NEPCore.compute_Mlincomb, Base.size
compute_Mder(::CustomNLEIGSNEP, λ::Number) = compute_Mder(sumnep, λ)
compute_Mlincomb(::CustomNLEIGSNEP, λ::Number, x) = compute_Mlincomb(sumnep, λ, x)
size(::CustomNLEIGSNEP, _) = n

@testset "NLEIGS: Custom NEP type" begin
    custom_nep = CustomNLEIGSNEP(n)
    @time X, lambda = nleigs(custom_nep, Sigma, maxit=10, v=ones(n), blksize=5)
    nleigs_verify_lambdas(4, custom_nep, X, lambda)
end
