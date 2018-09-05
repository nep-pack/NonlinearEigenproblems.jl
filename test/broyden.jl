using NonlinearEigenproblems.NEPSolver
using NonlinearEigenproblems.Gallery
using Test
using LinearAlgebra

@testset "broyden" begin
    dep=nep_gallery("dep0");
    S,V=broyden(dep)
    # Broyden returns a Schur factorization so check with MM
    @test opnorm(compute_MM(dep,S,V))<eps()*1000
    # test addconj
    S,V=broyden(dep,addconj=true,pmax=5)
    # Test by computing the eigenpairs
    λv,X = eigen(S)
    V_eigvecs=V*X;
    for k=1:size(λv,1)
        normalize!(V_eigvecs[:,k])
        @test norm(compute_Mlincomb(dep,λv[k],V_eigvecs[:,k]))<eps()*1000
    end
end
