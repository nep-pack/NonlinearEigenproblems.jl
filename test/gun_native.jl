# Intended to be run from nep-pack/ directory or nep-pack/test directory
if !isdefined(:global_modules_loaded)
    workspace()

    push!(LOAD_PATH, string(@__DIR__, "/../src"))

    using NEPCore
    using NEPTypes
    using LinSolvers
    using NEPSolver
    using Gallery
    using Base.Test
end

@testset "GUN (native)" begin
    nep = nep_gallery("nlevp_native_gun")
    n = size(nep, 1)
    tol = 1e-11
    ref_eigenvalue = 22345.116783765 + 0.644998598im # from NLEIGS

    @testset "Running algorithm" begin
        λ1,v1 = quasinewton(nep, λ = 150^2+1im, v = ones(n), displaylevel = 1, tol = tol, maxit = 500)

        v1 = v1 / norm(v1)

        @test norm(compute_Mlincomb(nep, λ1, v1)) < tol*100
        @test norm(compute_Mder(nep, λ1) * v1) < tol*100

        @test norm(λ1 - ref_eigenvalue) < tol*100
    end

    @testset "Compute derivatives" begin
        λ = 150^2+2im;
        V = randn(n, 2);
        z1 = compute_Mlincomb(nep, λ, V[:,1], [1.0], 1)

        # Compare with divided difference
        ee = 1e-4
        v = V[:,1];
        z2 = (compute_Mlincomb(nep, λ+ee, v) - compute_Mlincomb(nep, λ-ee, v)) / (2*ee);

        @test norm(z2-z1) < (ee^2)*1000
    end
end
