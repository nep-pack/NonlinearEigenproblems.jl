#  Tests for the projected NEPs
workspace()
push!(LOAD_PATH, pwd()*"/src")	

using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinSolvers

#using Winston # For plotting

using Base.Test


@testset "Rayleigh Function Iteration" begin
    println("Running two-sided RFI on random dep")
    nep=nep_gallery("dep0")
    nept=DEP([nep.A[1]',nep.A[2]'],copy(nep.tauv))

    n=size(nep,1);
    λ=NaN;
    x=NaN
    y=NaN;
    try
        λ,x,y =rfi(nep,nept,displaylevel=1,
                   v=ones(n), u=ones(n),
                   tolerance=1e-15);
    catch e
        # Only catch NoConvergence
        isa(e, NoConvergenceException) || rethrow(e)
        println("No convergence because:"*e.msg)
        # still access the approximations
        λ=e.λ
        x=e.v
    end
    println(λ)
    r1=compute_resnorm(nep,λ,x)
    println("Resnorm M:",r1)
    @test r1 < eps()*100

    r2=compute_resnorm(nept,λ,y)
    println("Resnorm M':",r2)
    @test r2 < eps()*100

    # Test RFIb

    λb,xb,yb =rfi_b(nep,nept,displaylevel=1,
                    v=v0, u=u0,λ=λ+0.01,tolerance=1e-15);
    @test λ ≈ λb

    
    println("Testing formula for derivative (with left and right eigvecs")
    tau=1;
    Mλ=-eye(n)-tau*nep.A[2]*exp(-tau*λ)
    Mtau= -λ*nep.A[2]*exp(-tau*λ);


    
    # Formula for derivative
    λp=-(y'*(Mtau)*x) / (y'* Mλ*x)

    δ=0.0001;



    nep.tauv[2]=nep.tauv[2]+δ
    nept.tauv[2]=nept.tauv[2]+δ

    λδ,x,y =rfi(nep,nept,displaylevel=1,
                v=ones(n), u=ones(n));

    λp_approx=(λδ-λ)/δ;

    @test abs(λp-λp_approx)< (δ*10)


    
end
