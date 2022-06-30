using NonlinearEigenproblems
using Test
using LinearAlgebra
using SparseArrays

@testset "ChebPEP" begin
    nep=nep_gallery("dep0",5);
    nep=DEP(nep.A./10);

    n=size(nep,1);
    chebpep=ChebPEP(nep,13,-1,1);

    # If we k is odd and the interval is symmetric, 0 is an
    # interpolation point
    s=0;
    @test norm(compute_Mder(nep,s)-compute_Mder(chebpep,s))<eps()*100

    # Check that we get a valid solution if we solve chebpep
    v0=ones(n);
    (λ2,V2)=iar(chebpep,neigs=2,errmeasure=ResidualErrmeasure(chebpep),v=v0)
    @test norm(compute_Mlincomb(nep,λ2[1],V2[:,1])) < eps()*5000

    # polyeig with chebyshev basis (non-standard interval)
    chebpep=ChebPEP(nep,16,-0.9,1.5);
    (λ,V)=polyeig(chebpep)
    j=argmin(abs.(λ));
    # Verify against the original problem
    @test norm(compute_Mlincomb(nep,λ[j],normalize(V[:,j]))) < eps()*1000
end
