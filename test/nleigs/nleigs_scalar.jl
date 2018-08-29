# Solves the scalar, purely nonlinear, eigenvalue problem
# A(λ) = 0.2*sqrt(λ) - 0.6*sin(2*λ) with both a polynomial
# and fully rational approach

# Intended to be run from nep-pack/ directory or nep-pack/test directory
if !@isdefined global_modules_loaded
    push!(LOAD_PATH, string(@__DIR__, "/../../src"))

    using NEPCore
    using NEPTypes
    using LinSolvers
    using NEPSolver
    using Gallery
    using IterativeSolvers
    using Test
end

include("nleigs_test_utils.jl")

as_matrix(x::Number) = (M = Matrix{eltype(x)}(1,1); M[1] = x; M)
n = 1
C = [as_matrix(0.2), as_matrix(-0.6)]
f = [λ -> sqrtm(λ), λ -> sin.(2*λ)]
nep = SPMF_NEP([C[1], C[2]], [f[1], f[2]])

Σ = complex([0.01, 4])

@testset "NLEIGS: Scalar (polynomial)" begin
    @time lambda, X = nleigs(nep, Σ, displaylevel=1, maxit=100, v=ones(n).+0im, leja=2, isfunm=false)

    # single eigenvalue converges
    nleigs_verify_lambdas(1, nep, X, lambda)
end

@testset "NLEIGS: Scalar (fully rational)" begin
    # set of poles candidates
    Ξ = -logspace(-6, 5, 10000)

    @time lambda, X = nleigs(nep, Σ, Ξ=Ξ, displaylevel=1, maxit=100, v=ones(n).+0im, leja=2, isfunm=false)

    # three eigenvalues converge
    nleigs_verify_lambdas(3, nep, X, lambda)
end
