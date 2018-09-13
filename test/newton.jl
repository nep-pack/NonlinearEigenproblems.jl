# Run tests on Newton methods & rfi methods

push!(LOAD_PATH, @__DIR__); using TestUtils

using NonlinearEigenproblems.NEPSolver
using NonlinearEigenproblems.NEPTypes
using NonlinearEigenproblems.Gallery
using Test
using LinearAlgebra

@testset "Newton iterations" begin
    nep=nep_gallery("dep0")

    # Construction to test the eigenvector extraction in the methods
    has_thrown_singexcep::Bool = false
    singexcep_errmeasure = function(λ,v)
        err = compute_resnorm(nep,λ,v)/norm(v)
        if (err < 1e-8)
            if has_thrown_singexcep
                err = 1
                global has_thrown_singexcep = false
            else
                global has_thrown_singexcep = true
                throw(SingularException(-1))
            end
        else
            global has_thrown_singexcep = false
        end
        return err
    end

    @bench @testset "Newton and AugNewton" begin
        println("Newton and AugNewton test")

        # newton and augnewton are equivalent, therefore I expect them
        # to generate identical results
        λ1,x1 = newton(nep,displaylevel=1,v=ones(size(nep,1)),λ=0,tol=eps()*10)
        λ2,x2 = augnewton(nep,displaylevel=1,v=ones(size(nep,1)),λ=0,tol=eps()*10)
        @test λ1 ≈ λ2
        @test x1 ≈ x2
        @test compute_resnorm(nep,λ1,x1) < eps()*100
        @test compute_resnorm(nep,λ2,x2) < eps()*100

        println("   eigenvector extraction")
        λ1,x1 = newton(nep,displaylevel=1,v=ones(size(nep,1)),λ=0,tol=eps()*10, errmeasure=singexcep_errmeasure)
        λ2,x2 = augnewton(nep,displaylevel=1,v=ones(size(nep,1)),λ=0,tol=eps()*10, errmeasure=singexcep_errmeasure)
        @test compute_resnorm(nep,λ1,x1) < 1e-8*100
        @test compute_resnorm(nep,λ2,x2) < 1e-8*100
    end

    @bench @testset "QuasiNewton" begin
    println("QuasiNewton  test")

        λ,x = quasinewton(nep,displaylevel=1,v=ones(size(nep,1)),λ=0,tol=1e-11)
        @test compute_resnorm(nep,λ,x) < 1e-11*100

        println("   eigenvector extraction")
        λ,x = quasinewton(nep,displaylevel=1,v=ones(size(nep,1)),λ=0,errmeasure=singexcep_errmeasure)
        @test compute_resnorm(nep,λ,x) < 1e-8*100

    end

    @bench @testset "Newton QR" begin

        println("Newton QR test")
        #Run derivative test for left and right eigenvectors returned by newtonqr
        λ3,x3,y3 =  newtonqr(nep, λ=0, v=ones(size(nep,1)), displaylevel=1,tol=eps()*10)
        @test compute_resnorm(nep,λ3,x3) < eps()*100

        println("   eigenvector extraction")
        λ3,x3,y3 =  newtonqr(nep, λ=0, v=ones(size(nep,1)), displaylevel=1,tol=eps()*10, errmeasure=singexcep_errmeasure)
        @test compute_resnorm(nep,λ3,x3) < 1e-8*100

        println("  Testing formula for derivative")
        tau = 1
        Mλ = -I - tau * nep.A[2] * exp(-tau*λ3)
        Mtau = -λ3*nep.A[2]*exp(-tau*λ3)

        # Formula for derivative
        λp=-(y3'*(Mtau)*x3) / (y3'* Mλ*x3)

        δ=0.0001
        nepp = deepcopy(nep)
        nepp.tauv[2] = nep.tauv[2]+δ
        λ3δ,x3,y3 = newtonqr(nepp, λ=0, v=ones(size(nep,1)), displaylevel=1,tol=eps()*10)
        λp_approx = (λ3δ-λ3)/δ
        @test abs(λp-λp_approx)< (δ*10)
    end

    @bench @testset "Resinv" begin
        println("resinv + periodicdde")
        nep1=nep_gallery("periodicdde", name="mathieu")
        # basic functionality of resinv. Start close to solution to speed up unit test
        λ4,x4 =  resinv(nep1, λ=-0.2447, v=[0.970208+0.0im, -0.242272+0.0im],
                        displaylevel=1,tol=eps()*10)
        r4=compute_resnorm(nep1,λ4,x4)

        @test r4 < eps()*100
    end


    @bench @testset "Rayleigh Function Iteration" begin
        println("rfi")
        nept = DEP([copy(nep.A[1]'), copy(nep.A[2]')],copy(nep.tauv))

        n=size(nep,1)
        u0=ones(n)
        v0=ones(n)


        λ,x,y =rfi(nep, nept, displaylevel=1, v=ones(n), u=ones(n), tol=1e-15)
        @test compute_resnorm(nep,λ,x) < eps()*100
        @test compute_resnorm(nept,λ,y) < eps()*100

        println("   eigenvector extraction")
        λ,x,y =rfi(nep, nept, displaylevel=1, v=ones(n), u=ones(n), tol=1e-15, errmeasure=singexcep_errmeasure)
        @test compute_resnorm(nep,λ,x) < 1e-8*100
        @test compute_resnorm(nept,λ,y) < 1e-8*100


        # Test RFIb

        println("rfi_b")
        λb,xb,yb =rfi_b(nep, nept, displaylevel=1, v=v0, u=u0, λ=λ+0.01, tol=1e-15)
        @test λ ≈ λb

        println("   eigenvector extraction")
        λb,xb,yb =rfi_b(nep, nept,displaylevel=1, v=v0, u=u0, λ=λ+0.01, tol=1e-15, errmeasure=singexcep_errmeasure)
        @test compute_resnorm(nep,λb,xb) < 1e-8*100
        @test compute_resnorm(nept,λb,yb) < 1e-8*100


        println("  Testing formula for derivative (with left and right eigvecs)")
        tau = 1
        Mλ = -I - tau*nep.A[2] * exp(-tau*λ)
        Mtau = -λ*nep.A[2]*exp(-tau*λ)

        # Formula for derivative
        λp=-(y'*(Mtau)*x) / (y'* Mλ*x)

        δ=0.0001
        nepp = deepcopy(nep)
        nepp.tauv[2] = nep.tauv[2]+δ
        neptp = deepcopy(nept)
        neptp.tauv[2] = nept.tauv[2]+δ
        λδ,x,y = rfi(nepp,neptp,displaylevel=1, v=ones(n), u=ones(n))
        λp_approx=(λδ-λ)/δ
        @test abs(λp-λp_approx)< (δ*10)



    end


    @bench @testset "implicitdet" begin
        println("Implicitdet test")
        nepd=nep_gallery("periodicdde", name="mathieu")
        λ,v=implicitdet(nepd, v=ones(size(nepd,1)), displaylevel=1)
        @test norm(compute_Mder(nepd,λ)*v) < eps()*100
    end
end
