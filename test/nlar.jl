# Unit test for the Nonlinear Arnoldi method (in src/method_nlar.jl)
# The "Gun" problem form gun_native.jl

# Intended to be run from nep-pack/ directory or nep-pack/test directory
push!(LOAD_PATH, string(@__DIR__, "/../src"))

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver
using Gallery
using LinearAlgebra
using Random
using Test

@testset "Nonlinear Arnoldi" begin
	TOL = 1e-10;

    println("\nTesting gun problem")
    Random.seed!(0)

    nep = nep_gallery("nlevp_native_gun")
    n = size(nep,1);

    #The eigenvalues of the gun problem lie within a semi-circle centred approximately around the "shift" with an approximate radius "scale"
    shift = 250^2;
    scale = 330^2-220^2;

    #Run Quasi-Newton for initial guess obtained from the knowledge of the eigenvalue distribution
    λ_ref,v_ref = quasinewton(nep, λ = shift+scale*(-0.131403+0.00759532im), v = ones(n), displaylevel = 2, tol = TOL/50, maxit = 500)
    println("Eigenvalue computed by quasi newton: ",λ_ref)

    #Shift and scale the NEP(mainly to avoid round-off errors because of the large entries in the coefficient matrices)
    nep1_spmf = SPMF_NEP(get_Av(nep),get_fv(nep))
    nep1 = shift_and_scale(nep1_spmf,shift=shift,scale=scale);

    #Run NLAR on the shifted and scaled NEP (nev set to 1 to save time. Method works for the case of finding multiple eigenvalues as well)
    λ,u = nlar(nep1, tol=TOL, λ0 = 0,maxit=500, nev = 4, R=0.01,mm =1,displaylevel=1,v0=ones(n),eigval_sorter=residual_eigval_sorter,inner_solver_method=NEPSolver.IARInnerSolver,qrfact_orth=false);
    λ_shifted = λ[1];v=u[:,1];

    #Rescale and shift back to the eigenvalue of the original problem
    λ_orig = shift+scale*λ_shifted;

    #Test residual and distance from reference eigenvalue
    @test norm(compute_Mlincomb(nep, λ_orig, v)) < TOL*100
    @test norm(compute_Mder(nep, λ_orig) * v) < TOL*100
    @test norm(λ_ref - λ_orig) < TOL*100
    ###########################################################################################################3

    #Testing PEP
    println("\nTesting a PEP")
    pepnep = nep_gallery("pep0",60)
    Dc,Vc = polyeig(pepnep,DefaultEigSolver)
    c = sortperm(abs.(Dc))
    println(" 6 smallest eigenvalues found by polyeig() according to the absolute values: \n ", Dc[c[1:6]])

    λ,u = nlar(pepnep, tol=TOL, maxit=60, nev = 3, R=0.001,mm=1,displaylevel=1, λ0=Dc[4]+1e-2,v0=ones(size(pepnep,1)), eigval_sorter=residual_eigval_sorter,inner_solver_method=NEPSolver.IARInnerSolver,qrfact_orth=false);
    println(" Smallest eigevalues found: \n λ: ",λ)

    # Test residuals
    @test norm(compute_Mlincomb(pepnep,λ[1],u[:,1])) < TOL
    @test norm(compute_Mlincomb(pepnep,λ[2],u[:,2])) < TOL
    @test norm(compute_Mlincomb(pepnep,λ[3],u[:,3])) < TOL

end
