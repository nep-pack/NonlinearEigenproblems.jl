# Run tests on Newton methods & rfi methods

# Intended to be run from nep-pack/ directory or nep-pack/test directory
if !isdefined(:global_running_all_tests) || global_running_all_tests != true
    workspace()

    push!(LOAD_PATH, string(@__DIR__, "/../src"))

    using NEPCore
    using NEPTypes
    using LinSolvers
    using NEPSolver
    using Gallery
    using Base.Test
end


nep=nep_gallery("dep0")
@testset "Newton iterations" begin

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

#Run derivative test for left and right eigenvectors returned by newtonqr
λ3,x3,y3 =  newtonqr(nep, λ=0, v=ones(size(nep,1)), displaylevel=0,tol=eps()*10);

#println("\nTesting formula for derivative (with left and right eigvecs)\n")
tau=1;
n=size(nep,1);
Mλ=-eye(n)-tau*nep.A[2]*exp(-tau*λ3)
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
# basic functionality of resinv
λ4,x4 =  resinv(nep, λ=0, v=ones(size(nep,1)), displaylevel=0,tol=eps()*10);
r4=compute_resnorm(nep,λ4,x4)

@test r4 < eps()*100
end

end


@testset "Rayleigh Function Iteration" begin
    println("Running two-sided RFI on random dep")
    nep=nep_gallery("dep0")
    nept=DEP([nep.A[1]',nep.A[2]'],copy(nep.tauv))

    n=size(nep,1);
    u0=ones(n);
    v0=ones(n);

    λ=NaN;
    x=NaN
    y=NaN;
    try
        λ,x,y =rfi(nep,nept,displaylevel=1,
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
    r1=compute_resnorm(nep,λ,x)
    println("Resnorm M:",r1)
    @test r1 < eps()*100

    r2=compute_resnorm(nept,λ,y)
    println("Resnorm M':",r2)
    @test r2 < eps()*100

    # Test RFIb

    λb,xb,yb =rfi_b(nep,nept,displaylevel=1,
                    v=v0, u=u0,λ=λ+0.01,tol=1e-15);
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
