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
        λ,W=ilan(nep,σ=0,γ=1;neigs=Inf,logger=0,maxit=30,tol=eps()*100,check_error_every=Inf,v=v0,errmeasure=ResidualErrmeasure(nep),inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=50))
        verify_lambdas(7, nep, λ, W, eps()*100)
    end

    @testset "Compute eigenpairs via Ritz extraction" begin
        nep=nep_gallery("dep_symm_double",10);
        v0=ones(size(nep,1));
        λ,W=ilan(nep;v=v0,tol=1e-5,neigs=3);
        verify_lambdas(3, nep, λ, W, 1e-5)
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

    @testset "Different format" begin
        n=100
        Random.seed!(1) # reset the random seed
        K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
        A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
        A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
        nep1=DEP([A1,A2],[0,1])
        nep2=SPMF_NEP([one(A1),A1,A2],[S -> -S,S -> one(S),S -> exp(-S)]);

        v0=rand(n)
        λ,W,err,V,H,ω,HH=ilan(nep1,σ=0,γ=1;neigs=Inf,logger=0,maxit=10,tol=eps()*100,check_error_every=Inf,v=v0)

        λ2,W2,err2,V2,H2,ω2,HH2,=ilan(nep2,σ=0,γ=1;neigs=Inf,logger=0,maxit=10,tol=eps()*100,check_error_every=Inf,v=v0)
        #@test opnorm(W - W2) < 1e-6
        @test norm(λ - λ)<1e-6
        @test norm(W - W2)<1e-6
        @test norm(V - V2)<1e-6
        @test norm(H - H2)<1e-6
        @test norm(ω - ω2)<1e-6
        @test norm(HH - HH2)<1e-6
        @test opnorm(V'*V - I) < 1e-6
        @test opnorm(V2'*V2 - I) < 1e-6



    end
end
