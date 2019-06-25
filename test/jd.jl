using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test
using LinearAlgebra

@testset "Jacobi–Davidson" begin

@test_throws ErrorException jd()


@bench @testset "Betcke-Voss" begin
@info "Test: Betcke-Voss"


@info "Testing a PEP"
nep = nep_gallery("pep0",60)
TOL = 1e-11;
λ,u = jd_betcke(nep, tol=TOL, maxit=55, neigs = 3, logger=displaylevel, v=ones(size(nep,1)),errmeasure=ResidualErrmeasure)
@info " Smallest eigenvalue found: $λ"
Dc,Vc = polyeig(nep,DefaultEigSolver)
c = sortperm(abs.(Dc))
@info " 6 smallest eigenvalues according to the absolute values: $(Dc[c[1:6]])"
verify_lambdas(3, nep, λ, u, TOL)


@info "Testing SG as inner solver"
nep = nep_gallery("real_quadratic")
nep = SPMF_NEP(get_Av(nep), get_fv(nep))
TOL = 1e-10;
λ,u=jd_betcke(Float64, nep, tol=TOL, maxit=4, logger = displaylevel, projtype = :Galerkin, inner_solver_method = SGIterInnerSolver, v=ones(size(nep,1)))
verify_lambdas(1, nep, λ, u, TOL)


@info "Testing IAR Cheb as projected solver"
nep = nep_gallery("dep0_sparse",40)
TOL = 1e-10;
λ,u = jd_betcke(ComplexF64, nep, tol=TOL, maxit=30, logger = displaylevel, inner_solver_method = IARChebInnerSolver, v=ones(size(nep,1)))
verify_lambdas(1, nep, λ, u, TOL)
λ0 = λ[1] # store these for the next test
v0 = vec(u)

@info "Testing convergence before starting"
λ,u=jd_betcke(nep, tol=TOL, maxit=25, neigs=1, logger=displaylevel, λ=λ0, v=v0)
verify_lambdas(1, nep, λ, u, TOL)


@info "Testing errors thrown"
nep = nep_gallery("pep0",4)
# Throw error if iterating more than the size of the NEP
@test_throws ErrorException λ,u=jd_betcke(nep, tol=TOL, maxit=60, logger = displaylevel, v=ones(size(nep,1)))
# SG requires Galerkin projection type to keep Hermitian
@test_throws ErrorException λ,u=jd_betcke(Float64, nep, tol=TOL, maxit=4, projtype = :PetrovGalerkin, inner_solver_method = SGIterInnerSolver, v=ones(size(nep,1)))
# An undefined projection type
@test_throws ErrorException λ,u=jd_betcke(nep, tol=TOL, maxit=4, projtype = :MYNOTDEFINED, v=ones(size(nep,1)))
# Too many required eigenvalues, will not converge and hence throw an exception
@test_throws NEPCore.NoConvergenceException λ,u=jd_betcke(nep, tol=TOL, maxit=4, neigs=1000, v=ones(size(nep,1)))

end



@bench @testset "Effenberger" begin
@info "Test: Effenberger"

TOL = 1e-10
nep = nep_gallery("pep0",250)
λ, u = jd_effenberger(nep, neigs=5, logger=displaylevel, tol=TOL, maxit=80, λ=0.82+0.9im, v=ones(ComplexF64,size(nep,1)))
verify_lambdas(5, nep, λ, u, TOL)

TOL = 1e-10
nep = nep_gallery("dep0",60)
λ, u = jd_effenberger(nep, neigs=3, logger=displaylevel, tol=TOL, maxit=55, λ=0.6+0im, v=ones(ComplexF64,size(nep,1)))#, inner_solver_method = IARChebInnerSolver)
verify_lambdas(3, nep, λ, u, TOL)

@info "Testing convergence before starting"
λ,u=jd_effenberger(nep, neigs=1, logger=displaylevel, tol=TOL, maxit=55, λ=λ[1], v=vec(u[:,1]))
verify_lambdas(1, nep, λ, u, TOL)


@info "Testing errors thrown"
nep = nep_gallery("pep0",50)
# Throw error if iterating more than the size of the NEP
@test_throws ErrorException λ, u = jd_effenberger(nep, tol=TOL, maxit=(size(nep,1)+1), v=ones(size(nep,1)))
# SG not possible with Effenberger
@test_throws ErrorException λ, u = jd_effenberger(nep, tol=TOL, maxit=40, inner_solver_method = SGIterInnerSolver, v=ones(size(nep,1)))
# Too many required eigenvalues, will not converge and hence throw an exception
@test_throws NEPCore.NoConvergenceException λ, u = jd_effenberger(nep, neigs=1000, tol=TOL, maxit=20, v=ones(size(nep,1)))

end

end
