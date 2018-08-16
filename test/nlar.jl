#Unit test for the Nonlinear Arnoldi method (in src/method_nlar.jl)

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

@testset "Nonlinear Arnoldi" begin
	TOL = 1e-10;

    println("\nTesting gun problem")
    nep = nep_gallery("nlevp_native_gun")
   	λ_ref = 22345.116783765 + 0.644998598im 			#from NLEIGS
   	λ0 = 150^2+1im;										#Initial guesses			
   	v0=compute_Mder(nep,λ0)\ones(size(nep,1))			#Perform one step of residual inverse iteration
	normalize!(v0)
    λ,u = nlar(nep, tol=TOL, λ0 = 150^2+1im,maxit=500, nev = 2, R=0.01,mm =1,displaylevel=1,v0=ones(size(nep,1)),inner_solver_method=NEPSolver.IARInnerSolver);
    λ1 = λ[1];v1=u[:,1];
    #Test distance from reference eigenvalue
    @test norm(compute_Mlincomb(nep, λ1, v1)) < tol*100
    @test norm(compute_Mder(nep, λ1) * v1) < tol*100
    @test norm(λ1 - λ_ref) < tol*100

    println("\nTesting a PEP")
	nep = nep_gallery("pep0",60)

	λ,u = nlar(nep, tol=TOL, maxit=60, nev = 3, R=0.1,mm =1,displaylevel=1,v=ones(size(nep,1)))
	println(" Smallest eigevalues found: \n λ: ",λ)
	Dc,Vc = polyeig(nep,DefaultEigSolver)
	c = sortperm(abs.(Dc))
	println(" 6 smallest eigenvalues according to the absolute values: \n ", Dc[c[1:6]])

	# Test residuals
	@test norm(compute_Mlincomb(nep,λ[1],u[:,1])) < TOL
	@test norm(compute_Mlincomb(nep,λ[2],u[:,2])) < TOL
	@test norm(compute_Mlincomb(nep,λ[3],u[:,3])) < TOL
end