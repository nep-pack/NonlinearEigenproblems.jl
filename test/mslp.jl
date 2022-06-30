using NonlinearEigenproblems
using Test
using LinearAlgebra
using SparseArrays

@testset "MSLP" begin

    target=0;
    nev=3;
    nep=nep_gallery("beam",5);

    # Full matrix test
    λ1,v1=mslp(nep)
    @test norm(compute_Mlincomb(nep,λ1,v1))<1e-10

    Av=[sparse(nep.A[1]),sparse(nep.A[2])];
    nep2 = DEP(nep.n, Av, nep.tauv)
    # Sparse matrix test
    λ2,v2 = mslp(nep2)
    @test norm(compute_Mlincomb(nep2,λ2,v2)) < 1e-10

    @test λ1≈λ2 # They should be exactly the same

    @info "mslp noconv + double"
    @bench @testset  "mslp + double" begin
        nep3 = nep_gallery("dep_double");
        @test_throws NoConvergenceException mslp(nep3, λ=9im, maxit=10)
        λ,v = mslp(nep3, λ=9im, maxit=100,errmeasure=ResidualErrmeasure(nep3))
        @test norm(compute_Mlincomb(nep3, λ, v)) < eps()*100
    end
end
