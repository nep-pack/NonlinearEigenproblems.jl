# Solves the scalar, purely nonlinear, eigenvalue problem
# A(λ) = 0.2*sqrt(λ) - 0.6*sin(2*λ) with both a polynomial
# and fully rational approach

#using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test

function nleigs_scalar()
    as_matrix(x::Number) = (M = Matrix{eltype(x)}(undef,1,1); M[1] = x; M)
    n = 1
    C = [as_matrix(0.2), as_matrix(-0.6)]
    f = [λ -> sqrt(λ), λ -> sin(2*λ)]
    nep = SPMF_NEP([C[1], C[2]], [f[1], f[2]])

    Σ = complex([0.01, 4])

    @bench @testset "Polynomial" begin
        lambda, X = nleigs(nep, Σ, logger=displaylevel, maxit=100, v=ones(n).+0im, leja=2, isfunm=false)

        # single eigenvalue converges
        verify_lambdas(1, nep, lambda, X)
    end

    @bench @testset "Fully rational" begin
        # set of poles candidates
        Ξ = -10 .^ range(-6, stop = 5, length = 10000)

        lambda, X = nleigs(nep, Σ, Ξ=Ξ, logger=displaylevel, maxit=100, v=ones(n).+0im, leja=2, isfunm=false)

        # three eigenvalues converge
        verify_lambdas(3, nep, lambda, X)
    end
end

@testset "NLEIGS: Scalar" begin
    nleigs_scalar()
end
