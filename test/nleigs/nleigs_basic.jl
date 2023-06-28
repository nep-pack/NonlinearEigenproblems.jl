# Solves a few basic eigenvalue problems to test various aspects of NLEIGS

#using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using NonlinearEigenproblems.RKHelper
using Test
using LinearAlgebra

function nleigs_basic()
    n = 2
    B = Vector{Matrix{Float64}}([[1 3; 5 6], [3 4; 6 6], [1 0; 0 1]])
    pep = PEP(B)

    Σ = [-10.0-2im, 10-2im, 10+2im, -10+2im]

    @bench @testset "Polynomial only" begin
        lambda, X = nleigs(pep, Σ, maxit=10, v=ones(n).+0im, blksize=5)
        verify_lambdas(4, pep, lambda, X)
    end

    @bench @testset "Non-convergent linearization" begin
        @test_logs (:warn, r".*Linearization not converged.*") match_mode = :any begin
            lambda, X = nleigs(pep, Σ, maxit=10, v=ones(n).+0im, maxdgr=5, blksize=5)
            verify_lambdas(4, pep, lambda, X)
        end
    end

    @bench @testset "Non-convergent linearization (static)" begin
        @test_logs (:warn, r".*Linearization not converged.*") match_mode = :any begin
            lambda, X = nleigs(pep, Σ, maxit=10, v=ones(n).+0im, maxdgr=5, blksize=5, static=true)
            verify_lambdas(4, pep, lambda, X)
        end
    end

    @bench @testset "Non-convergent linearization (return_details)" begin
        @test_logs (:warn, r".*Linearization not converged.*") match_mode = :any begin
            lambda, X, _ = nleigs(pep, Σ, maxit=5, v=ones(n).+0im, blksize=5, return_details=true)
            verify_lambdas(0, pep, lambda, X)
        end
    end

    @bench @testset "Complex-valued matrices" begin
        complex_B = map(X -> X + im*I, B)
        complex_pep = PEP(complex_B)
        lambda, X, _ = nleigs(complex_pep, Σ, maxit=10, v=ones(n).+0im, blksize=5, return_details=true)
        verify_lambdas(3, complex_pep, lambda, X)
    end

    @bench @testset "Complex-valued start vector" begin
        lambda, X, _ = nleigs(pep, Σ, maxit=10, v=ones(n) * (1+0.1im), blksize=5, return_details=true)
        verify_lambdas(4, pep, lambda, X)
    end

    @bench @testset "return_details" begin
        lambda, X, res, details = nleigs(pep, Σ, maxit=10, v=ones(n).+0im, blksize=5, return_details=true)
        verify_lambdas(4, pep, lambda, X)

        lam = details.Lam[:,end]
        res = details.Res[:,end]
        conv = map(i -> res[i] < 1e-12 && inpolygon(real(lam[i]), imag(lam[i]), real(Σ), imag(Σ)), 1:size(lam, 1))
        lamconv = lam[conv, end]

        # test that eigenvalues in the info are the same as those returned by nleigs
        @test length(lamconv) == 4
        @test length(union(lambda, lamconv)) == 4
    end
end

@testset "NLEIGS: Basic functionality" begin
    nleigs_basic()
end
