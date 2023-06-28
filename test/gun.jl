# Run tests on Beyns contour integral method

using NonlinearEigenproblems
using Test
using LinearAlgebra

push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
using GalleryNLEVP

@testset "GUN (NLEVP interface)" begin

    nep_org = nep_gallery(NLEVP_NEP,"gun")
    nep1 = nep_gallery("nlevp_native_gun")
    n=size(nep1,1);
    tol=1e-11;
    @bench @testset "Running alg" begin
        λ1,v1=quasinewton(nep1,λ=150^2+1im,v=ones(n),logger=displaylevel,tol=tol,maxit=500);
        λ1 = λ1[1]
        v1 = vec(v1)
        normalize!(v1)

        @test norm(compute_Mlincomb(nep1,λ1,v1))<tol*100
        @test norm(compute_Mder(nep1,λ1)*v1)<tol*100

        @test norm(compute_Mlincomb(nep_org,λ1,v1))<tol*100
        @test norm(compute_Mder(nep_org,λ1)*v1)<tol*100

    end
    @bench @testset "Compute derivatives" begin
        # Test compute_Mlincomb:
        λ=150^2+2im;
        V=randn(n,2);    v=V[:,1];
        z1=compute_Mlincomb(nep1,λ,V[:,1], [1.0], 1)
        z2=compute_Mlincomb(nep_org,λ,V[:,1], [1.0], 1)
        # Compare with divided difference
        ee=1e-4
        z3=(compute_Mlincomb(nep_org,λ+ee,v)-compute_Mlincomb(nep_org,λ-ee,v))/(2*ee);

        @test norm(z3-z1)<(ee^2)*1000  # Compare native and FD
        @test norm(z3-z2)<(ee^2)*1000 # Compare MATLAB and FD. This fails on NLEVP 3.0 due to a bug in gun.m in NLEVP
        @test norm(z1-z2)<sqrt(eps()) # Compare MATLAB and native.  This fails on NLEVP 3.0 due to a bug in gun.m in NLEVP
    end


    @bench @testset "Compare SPMF and Native in alg" begin
        # Check that two steps of quasinewton always gives the same result
        λ_org=0
        try
            quasinewton(nep_org,maxit=2,λ=150^2+1im,v=ones(n),logger=displaylevel)
        catch e
            λ_org=e.λ
        end


        λ1=0
        try
            quasinewton(nep1,maxit=2,λ=150^2+1im,v=ones(n),logger=displaylevel)
        catch e
            λ1=e.λ
        end

        @test abs(λ1-λ_org)<sqrt(eps())

    end

    @bench @testset "Compare MATLAB loaded vs MAT-loaded" begin
        nep2=nep_gallery("nlevp_native_gun");
        z=ones(n); λ=150^2;
        @test (norm(compute_Mlincomb(nep2,λ,z)-compute_Mlincomb(nep1,λ,z))/norm(compute_Mlincomb(nep1,λ,z)) < 1e-14)
    end

end
