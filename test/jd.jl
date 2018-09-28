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
λ,u = jd_betcke(nep, tol=TOL, maxit=55, Neig = 3, displaylevel=displaylevel, v=ones(size(nep,1)))
@info " Smallest eigenvalue found: $λ"
Dc,Vc = polyeig(nep,DefaultEigSolver)
c = sortperm(abs.(Dc))
@info " 6 smallest eigenvalues according to the absolute values: $(Dc[c[1:6]])"

# Test residuals
@test norm(compute_Mlincomb(nep,λ[1],u[:,1])) < TOL
@test norm(compute_Mlincomb(nep,λ[2],u[:,2])) < TOL
@test norm(compute_Mlincomb(nep,λ[3],u[:,3])) < TOL



@info "Testing SG as inner solver"
nep = nep_gallery("real_quadratic")
nep = SPMF_NEP(get_Av(nep), get_fv(nep))
TOL = 1e-10;
# Also test that a warning is issued
λ,u=jd_betcke(Float64, nep, tol=TOL, maxit=4, displaylevel = displaylevel, projtype = :Galerkin, inner_solver_method = NEPSolver.SGIterInnerSolver, v=ones(size(nep,1)))
λ = λ[1]
u = vec(u)
@info " Resnorm of computed solution: $(compute_resnorm(nep,λ,u))"
@info " Smallest eigenvalue found: $λ"

@test norm(compute_Mlincomb(nep,λ,u)) < TOL



@info "Testing IAR Cheb as projected solver"
nep = nep_gallery("dep0_sparse",40)
TOL = 1e-10;
λ,u = jd_betcke(ComplexF64, nep, tol=TOL, maxit=30, displaylevel = displaylevel, inner_solver_method = NEPSolver.IARChebInnerSolver, v=ones(size(nep,1)))
λ = λ[1]
u = vec(u)
@info " Resnorm of computed solution: $(compute_resnorm(nep,λ,u))"
@info " Smallest eigenvalue found: $λ"

@test norm(compute_Mlincomb(nep,λ,u)) < TOL



@info "Testing convergence before starting"
λ,u=jd_betcke(nep, tol=TOL, maxit=25, Neig=1, displaylevel=displaylevel, λ=λ, v=u)
λ = λ[1]
u = vec(u)
@test norm(compute_Mlincomb(nep,λ,u)) < TOL



@info "Testing errors thrown"
nep = nep_gallery("pep0",4)
# Throw error if iterating more than the size of the NEP
@test_throws ErrorException λ,u=jd_betcke(nep, tol=TOL, maxit=60, displaylevel = displaylevel, v=ones(size(nep,1)))
# SG requires Galerkin projection type to keep Hermitian
@test_throws ErrorException λ,u=jd_betcke(Float64, nep, tol=TOL, maxit=4, projtype = :PetrovGalerkin, inner_solver_method = NEPSolver.SGIterInnerSolver, v=ones(size(nep,1)))
# An undefined projection type
@test_throws ErrorException λ,u=jd_betcke(nep, tol=TOL, maxit=4, projtype = :MYNOTDEFINED, v=ones(size(nep,1)))
# Too many required eigenvalues, will not converge and hence throw an exception
@test_throws NEPCore.NoConvergenceException λ,u=jd_betcke(nep, tol=TOL, maxit=4, Neig=1000, v=ones(size(nep,1)))

end



@bench @testset "Effenberger" begin
@info "Test: Effenberger"

TOL = 1e-10
nep = nep_gallery("pep0",60)
λ, u = jd_effenberger(nep, Neig=3, displaylevel=displaylevel, tol=TOL, maxit=55, λ=0.82+0.9im, v=ones(ComplexF64,size(nep,1)))
@info " Eigenvalues found: $λ"
@test norm(compute_Mlincomb(nep,λ[1],u[:,1])) < TOL
@test norm(compute_Mlincomb(nep,λ[2],u[:,2])) < TOL
@test norm(compute_Mlincomb(nep,λ[3],u[:,3])) < TOL

TOL = 1e-10
nep = nep_gallery("dep0",60)
λ, u = jd_effenberger(nep, Neig=3, displaylevel=displaylevel, tol=TOL, maxit=55, λ=0.6+0im, v=ones(ComplexF64,size(nep,1)))#, inner_solver_method = NEPSolver.IARChebInnerSolver)
@info " Eigenvalues found: $λ"
@test norm(compute_Mlincomb(nep,λ[1],u[:,1])) < TOL
@test norm(compute_Mlincomb(nep,λ[2],u[:,2])) < TOL
@test norm(compute_Mlincomb(nep,λ[3],u[:,3])) < TOL

@info "Testing convergence before starting"
λ,u=jd_effenberger(nep, Neig=1, displaylevel=displaylevel, tol=TOL, maxit=55, λ=λ[1], v=vec(u[:,1]))
λ = λ[1]
u = vec(u)
@test norm(compute_Mlincomb(nep,λ,u)) < TOL


@info "Testing errors thrown"
nep = nep_gallery("pep0",50)
# Throw error if iterating more than the size of the NEP
@test_throws ErrorException λ, u = jd_effenberger(nep, tol=TOL, maxit=(size(nep,1)+1), v=ones(size(nep,1)))
# SG not possible with Effenberger
@test_throws ErrorException λ, u = jd_effenberger(nep, tol=TOL, maxit=40, inner_solver_method = NEPSolver.SGIterInnerSolver, v=ones(size(nep,1)))
# Too many required eigenvalues, will not converge and hence throw an exception
@test_throws NEPCore.NoConvergenceException λ, u = jd_effenberger(nep, Neig=1000, tol=TOL, maxit=20, v=ones(size(nep,1)))

end

end
