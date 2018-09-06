#  Tests for the Linear solvers

using Test
using LinearAlgebra

using LinSolversMATLAB

@testset "LinSolvers" begin
    A=sparse(randn(5,5));

    ## Eigs solvers (native and MATLAB)
    DD1=NativeEigSSolver(A)
    λ,V=eig_solve(DD1,nev=3)

    λ1=λ[1];
    v1=V[:,1];
    @test opnorm(A*v1-λ1*v1)/opnorm(v1)<eps()*100

    DD2=MatlabEigSSolver(A)
    λ,V=eig_solve(DD2,nev=3)


    λ1=λ[1];
    v1=V[:,1];
    @test opnorm(A*v1-λ1*v1)/opnorm(v1)<eps()*100


    ## Eigs solvers for GEP (only MATLAB since julia buggy)
    B=sparse(randn(5,5));


    DD2=MatlabEigSSolver(A,B)
    λ,V=eig_solve(DD2,nev=3)

    λ1=λ[1];
    v1=V[:,1];
    @test opnorm(A*v1-λ1*B*v1)/opnorm(v1)<eps()*100
end
