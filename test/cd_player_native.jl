using NonlinearEigenproblems
using Test
using LinearAlgebra


@bench @testset "cd player native" begin
    nep = nep_gallery("nlevp_native_cd_player")
    n = size(nep, 1)

    (Λ,V)=polyeig(nep)

    (λ,v)=quasinewton(nep,λ=Λ[1],v=V[:,1],logger=displaylevel, armijo_factor=0.5,armijo_max=10, tol=1e-10)
    verify_lambdas(1, nep, λ, v, 1e-10)

    (λ,v)=resinv(nep,λ=Λ[2],v=V[:,2],logger=displaylevel, armijo_factor=0.5,armijo_max=10, tol=1e-10)
    verify_lambdas(1, nep, λ, v, 1e-10)

    @testset "Errors thrown" begin
        @test_throws MethodError nep_gallery("nlevp_native_cd_player", 15)
        @test_throws MethodError nep_gallery("nlevp_native_cd_player", t=15)
    end

end
