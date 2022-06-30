# Unit test for the Nonlinear Arnoldi method (in src/method_nlar.jl)
# The "Gun" problem form gun_native.jl

using NonlinearEigenproblems
using Test
using LinearAlgebra
using Random

@bench @testset "Nonlinear Arnoldi" begin
	TOL = 1e-10;

    @info "Testing gun problem"
    Random.seed!(0)

    nep = nep_gallery("nlevp_native_gun")
    n = size(nep,1);

    #The eigenvalues of the gun problem lie within a semi-circle centred approximately around the "shift" with an approximate radius "scale"
    shift = 250^2;
    scale = 330^2-220^2;

    #Run Quasi-Newton for initial guess obtained from the knowledge of the eigenvalue distribution
    λ_ref,v_ref = quasinewton(nep, λ = shift+scale*(-0.131403+0.00759532im), v = ones(n), logger = displaylevel, tol = TOL/50, maxit = 500)
    @info "Eigenvalue computed by quasi newton: $λ_ref"

    #Shift and scale the NEP(mainly to avoid round-off errors because of the large entries in the coefficient matrices)
    nep1_spmf = SPMF_NEP(get_Av(nep),get_fv(nep))
    nep1 = shift_and_scale(nep1_spmf,shift=shift,scale=scale);

    #Run NLAR on the shifted and scaled NEP (neigs set to 1 to save time. Method works for the case of finding multiple eigenvalues as well)
    λ,u = nlar(nep1, tol=TOL, λ=0, maxit=100, neigs=2, R=0.01, logger=displaylevel, v=ones(n),
        eigval_sorter=residual_eigval_sorter, inner_solver_method=IARInnerSolver(),
        qrfact_orth=false, num_restart_ritz_vecs=8, max_subspace=150)

    λ_shifted = λ[1];v=u[:,1];

    #Rescale and shift back to the eigenvalue of the original problem
    λ_orig = shift+scale*λ_shifted;

    #Test residual and distance from reference eigenvalue
    @test norm(compute_Mlincomb(nep, λ_orig, v)) < sqrt(TOL)*50
    @test norm(compute_Mder(nep, λ_orig) * v) < sqrt(TOL)*50
    @test norm(λ_ref - λ_orig) < sqrt(TOL)*500
    ###########################################################################################################3

    #Testing PEP
    @info "Testing a PEP"
    pepnep = nep_gallery("pep0",60)
    Dc,Vc = polyeig(pepnep,DefaultEigSolver)
    c = sortperm(abs.(Dc))
    @info " 6 smallest eigenvalues found by polyeig() according to the absolute values: $(Dc[c[1:6]])"

    λ,u = nlar(pepnep, tol=TOL, maxit=60, neigs=3, R=0.001, logger=displaylevel,
        λ=Dc[4]+1e-2, v=ones(size(pepnep,1)), eigval_sorter=residual_eigval_sorter,
        inner_solver_method=IARInnerSolver(), qrfact_orth=false, errmeasure=ResidualErrmeasure(pepnep))

    verify_lambdas(3, pepnep, λ, u, TOL)


    @info "Testing errors thrown"
    nep = nep_gallery("pep0",4)
    @test_throws NEPCore.NoConvergenceException λ,u= nlar(nep, tol=1e-20, maxit=2, neigs=3)

end
