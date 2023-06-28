using NonlinearEigenproblems
using Test
using LinearAlgebra


@bench @testset "fiber native" begin
    nep = nep_gallery("nlevp_native_fiber")
    n = size(nep, 1)

    # An exact eigenvalue according (reported in NLEVP collection)
    sol_val= 7.139494306065948e-07;

    (λ,v)=quasinewton(nep,λ=7.14e-7,v=ones(n),
                      logger=displaylevel, armijo_factor=0.5,armijo_max=10)

    @test abs(λ-sol_val)<1e-10;

    # check that we "maintain" real arithmetic
    vv=real(v/v[1]);
    (λ1,v)=resinv(Float64,nep,λ=7.14e-7,v=vv,
                      logger=displaylevel)
    @test abs(λ-sol_val)<1e-10;

    @test eltype(v)==Float64

    @testset "Errors thrown" begin
        @test_throws MethodError nep_gallery("nlevp_native_fiber", 15)
        @test_throws MethodError nep_gallery("nlevp_native_fiber", t=15)
    end

end
