using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test
using IterativeSolvers
using LinearAlgebra
using Random
using SparseArrays
using Revise
include("../src/method_ilan.jl");


@testset "ILAN" begin

    n=1000
    Random.seed!(1) # reset the random seed
    K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
    A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
    A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
    A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';

    nep=DEP([A1,A2,A3],[0,1,0.8])
    v0=rand(n)
    V,H,ω,HH=ilan_benchmark(nep,σ=0,γ=1;Neig=10,displaylevel=1,maxit=10,tol=eps()*100,check_error_every=1,v=v0)
    V2,H2,ω2,HH2=ilan(nep,σ=0,γ=1;Neig=10,displaylevel=1,maxit=10,tol=eps()*100,check_error_every=1,v=v0)

    @test norm(V-V2)<sqrt(eps())
    @test norm(H-H2)<sqrt(eps())
    @test norm(ω-ω2)<sqrt(eps())
    @test norm(HH-HH2)<sqrt(eps())


end
