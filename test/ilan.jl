using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test
using IterativeSolvers
using LinearAlgebra
using Random
using SparseArrays


@testset "ILAN" begin

    @testset "Compute as many eigenpairs as possible (neigs=Inf)" begin
        n=100
        Random.seed!(1) # reset the random seed
        K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
        A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
        A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
        A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';
        nep=DEP([A1,A2,A3],[0,1,0.8])
        v0=rand(n)
        λ,W,V2,H2,ω2,HH2=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=30,tol=eps()*100,check_error_every=5,v=v0,errmeasure=ResidualErrmeasure(nep))
        verify_lambdas(3, nep, λ, W, eps()*100)
    end

    @testset "Errors thrown" begin
        n=100
        Random.seed!(1) # reset the random seed
        K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
        A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
        A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
        A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';
        nep=DEP([A1,A2,A3],[0,1,0.8])
        v0=rand(n)
        @test_throws NEPCore.NoConvergenceException λ,W,V2,H2,ω2,HH2=ilan(nep,σ=0,γ=1;neigs=2,logger=0,maxit=3,tol=eps()*100,check_error_every=Inf,v=v0,errmeasure=ResidualErrmeasure(nep))
    end

end
