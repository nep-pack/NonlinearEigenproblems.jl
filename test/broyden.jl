using NonlinearEigenproblems
using Test
using LinearAlgebra

@bench @testset "broyden" begin
    dep=nep_gallery("dep1");
    S,V=broyden(dep)
    # Broyden returns a Schur factorization so check with MM
    @test opnorm(compute_MM(dep,S,V))<sqrt(eps())
    # test addconj
    S,V=broyden(dep,addconj=true,pmax=5)
    # Test by computing the eigenpairs
    λv,X = eigen(S)
    V_eigvecs=V*X;
    verify_lambdas(4, dep, λv, V_eigvecs, sqrt(eps()))
end
