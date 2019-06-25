using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test
using IterativeSolvers
using LinearAlgebra
using Random
using SparseArrays


@testset "ILAN" begin

    @testset "DEP" begin
        n=1000
        Random.seed!(1) # reset the random seed
        K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
        A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
        A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
        A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';

        nep=DEP([A1,A2,A3],[0,1,0.8])
        v0=rand(n)
        V,H,ω,HH=ilan_benchmark(nep,σ=0,γ=1;neigs=10,displaylevel=0,maxit=10,tol=eps()*100,check_error_every=1,v=v0,errmeasure=ResidualErrmeasure)
        _,_,V2,H2,ω2,HH2=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=10,tol=eps()*100,check_error_every=1,v=v0,errmeasure=ResidualErrmeasure)

        # Disabled. This unit test requires rewriting (and ilan_benchmark can be removed)
        @test norm(V-V2)<sqrt(eps())
        @test norm(H-H2)<sqrt(eps())
        @test norm(ω-ω2)<sqrt(eps())
        @test norm(HH-HH2)<sqrt(eps())
    end

    @testset "SPMF and DEP" begin

        n=1000
        Random.seed!(1) # reset the random seed
        K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
        A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
        A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
        A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';

        nep=DEP([A1,A2,A3],[0,1,0.8])
        v0=rand(n)
        V,H,ω,HH=ilan_benchmark(nep,σ=0,γ=1;neigs=10,displaylevel=0,maxit=10,tol=eps()*100,check_error_every=Inf,v=v0,errmeasure=ResidualErrmeasure)

        f1= S -> -S; f2= S -> one(S); f3= S -> exp(-S); f4= S -> exp(-0.8*S)
        nep=SPMF_NEP([one(A1),A1,A2,A3],[f1,f2,f3,f4])
        _,_,V2,H2,ω2,HH2=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=10,tol=eps()*100,check_error_every=Inf,v=v0,errmeasure=ResidualErrmeasure)

        @test norm(V-V2)<sqrt(eps())
        @test norm(H-H2)<sqrt(eps())
        @test norm(ω-ω2)<sqrt(eps())
        @test norm(HH-HH2)<sqrt(eps())
    end

    @testset "Compute as many eigenpairs as possible (neigs=Inf)" begin
        n=100
        Random.seed!(1) # reset the random seed
        K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
        A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
        A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
        A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';
        nep=DEP([A1,A2,A3],[0,1,0.8])
        v0=rand(n)
        λ,W,V2,H2,ω2,HH2=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=30,tol=eps()*100,check_error_every=5,v=v0,errmeasure=ResidualErrmeasure)
        verify_lambdas(2, nep, λ, W, eps()*100)
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
        @test_throws NEPCore.NoConvergenceException λ,W,V2,H2,ω2,HH2=ilan(nep,σ=0,γ=1;neigs=2,logger=0,maxit=3,tol=eps()*100,check_error_every=Inf,v=v0,errmeasure=ResidualErrmeasure)
    end

end
