# Run tests on Newton methods & rfi methods

using NonlinearEigenproblems.NEPSolver
using NonlinearEigenproblems.NEPTypes
using NonlinearEigenproblems.Gallery
using Test
using LinearAlgebra

@testset "Newton iterations" begin
    nep=nep_gallery("dep0")

    @testset "Newton and AugNewton" begin
        # newton and augnewton are equivalent, therefore I expect them
        # to generate identical results
        λ1,x1 =newton(nep,displaylevel=0,v=ones(size(nep,1)),λ=0,tol=eps()*10);
        λ2,x2 =augnewton(nep,displaylevel=0,v=ones(size(nep,1)),λ=0,tol=eps()*10);

        @test λ1 ≈ λ2
        @test x1 ≈ x2
        r1=compute_resnorm(nep,λ1,x1)
        r2=compute_resnorm(nep,λ2,x2)

        @test r1 < eps()*100
        @test r2 < eps()*100
    end

    @testset "Newton QR" begin

        println("Newton QR test")
        #Run derivative test for left and right eigenvectors returned by newtonqr
        λ3,x3,y3 =  newtonqr(nep, λ=0, v=ones(size(nep,1)), displaylevel=0,tol=eps()*10);

        #println("\nTesting formula for derivative (with left and right eigvecs)\n")
        tau=1;
        Mλ = -I - tau * nep.A[2] * exp(-tau*λ3)
        Mtau= -λ3*nep.A[2]*exp(-tau*λ3);

        # Formula for derivative
        λp=-(y3'*(Mtau)*x3) / (y3'* Mλ*x3)

        δ=0.0001;

        nep.tauv[2]=nep.tauv[2]+δ

        λ3δ,x3,y3 =newtonqr(nep, λ=0, v=ones(size(nep,1)), displaylevel=0,tol=eps()*10);

        λp_approx=(λ3δ-λ3)/δ;

        @test abs(λp-λp_approx)< (δ*10)
    end

    @testset "Resinv" begin
        println("resinv + periodicdde")
        nep1=nep_gallery("periodicdde","rand0")
        # basic functionality of resinv. Start close to solution to speed up unit test
        λ4,x4 =  resinv(nep1, λ=-0.2447, v=[0.970208+0.0im, -0.242272+0.0im],
                        displaylevel=1,tol=eps()*10);
        r4=compute_resnorm(nep1,λ4,x4)

        @test r4 < eps()*100
    end


    @testset "Rayleigh Function Iteration" begin
        println("Running two-sided RFI on random dep")
        nepd=nep_gallery("dep0")
        nept=DEP([nepd.A[1]',nepd.A[2]'],copy(nepd.tauv))

        n=size(nepd,1);
        u0=ones(n);
        v0=ones(n);

        λ=NaN;
        x=NaN
        y=NaN;
        try
            λ,x,y =rfi(nepd,nept,displaylevel=1,
                       v=ones(n), u=ones(n),
                       tol=1e-15);
        catch e
            # Only catch NoConvergence
            isa(e, NoConvergenceException) || rethrow(e)
            println("No convergence because:"*e.msg)
            # still access the approximations
            λ=e.λ
            x=e.v
        end
        println(λ)
        r1=compute_resnorm(nepd,λ,x)
        println("Resnorm M:",r1)
        @test r1 < eps()*100

        r2=compute_resnorm(nept,λ,y)
        println("Resnorm M':",r2)
        @test r2 < eps()*100

        # Test RFIb

        λb,xb,yb =rfi_b(nepd,nept,displaylevel=1,
                        v=v0, u=u0,λ=λ+0.01,tol=1e-15);
        @test λ ≈ λb


        println("Testing formula for derivative (with left and right eigvecs")
        tau=1;
        Mλ = -I - tau*nepd.A[2] * exp(-tau*λ)
        Mtau= -λ*nepd.A[2]*exp(-tau*λ);


        # Formula for derivative
        λp=-(y'*(Mtau)*x) / (y'* Mλ*x)

        δ=0.0001;


        nepd.tauv[2]=nepd.tauv[2]+δ
        nept.tauv[2]=nept.tauv[2]+δ

        λδ,x,y =rfi(nepd,nept,displaylevel=1,
                    v=ones(n), u=ones(n));

        λp_approx=(λδ-λ)/δ;

        @test abs(λp-λp_approx)< (δ*10)



    end


    @testset "implicitdet" begin
        nepd=nep_gallery("periodicdde","mathieu")
        λ,v=implicitdet(nepd, v=ones(size(nepd,1)))
        @test norm(compute_Mder(nepd,λ)*v) < eps()*100
    end
end
