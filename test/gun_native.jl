using NonlinearEigenproblems
using Test
using LinearAlgebra

@testset "GUN (native)" begin
    nep = nep_gallery("nlevp_native_gun")
    n = size(nep, 1)
    tol = 1e-11
    ref_eigenvalue = 22345.116783765 + 0.644998598im # from NLEIGS

    @bench @testset "Running algorithm" begin
        λ1,v1 = quasinewton(nep, λ = 150^2+1im, v = ones(n), logger = displaylevel, tol = tol, maxit = 500)

        v1 = v1 / norm(v1)

        @test norm(compute_Mlincomb(nep, λ1, v1)) < sqrt(tol)
        @test norm(compute_Mder(nep, λ1) * v1) < sqrt(tol)

        @test norm(λ1 - ref_eigenvalue) < sqrt(tol)*100
    end

    @bench @testset "Compute derivatives" begin
        λ = 150^2+2im
        v = randn(n)
        z1 = compute_Mlincomb(nep, λ, v, [1.0], 1)

        # Compare with divided difference
        ee = 1e-4
        z2 = (compute_Mlincomb(nep, λ+ee, v) - compute_Mlincomb(nep, λ-ee, v)) / (2*ee)

        @test norm(z2-z1) < (ee^2)*1000
    end

    @testset "Errors thrown" begin
        @test_throws MethodError nep_gallery("nlevp_native_gun", 15)
        @test_throws MethodError nep_gallery("nlevp_native_gun", t=15)
    end
end
