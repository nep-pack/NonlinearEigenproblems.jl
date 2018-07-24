if !isdefined(:global_running_all_tests) || global_running_all_tests != true
    workspace()

    push!(LOAD_PATH, string(@__DIR__, "/../src"))

    using NEPCore
    using NEPTypes
    using LinSolvers
    using NEPSolver
    using Gallery
    using IterativeSolvers
    using Base.Test
end


@testset "Jacobi–Davidson" begin



println("\nTesting a PEP")
nep = nep_gallery("pep0",60)
TOL = 1e-10;
λ,u = jd(nep, tol=TOL, maxit=55, Neig = 3, displaylevel=1, v0=ones(size(nep,1)))
println(" Smallest eigevalues found: \n λ: ",λ)
Dc,Vc = polyeig(nep,DefaultEigSolver)
c = sortperm(abs.(Dc))
println(" 6 smallest eigenvalues according to the absolute values: \n ", Dc[c[1:6]])

# Test residuals
@test norm(compute_Mlincomb(nep,λ[1],u[:,1])) < TOL
@test norm(compute_Mlincomb(nep,λ[2],u[:,2])) < TOL
@test norm(compute_Mlincomb(nep,λ[3],u[:,3])) < TOL

# Test realtive errors in eig-values and take hight for complex conjuagte values
@test  abs(Dc[c[1]]-λ[1])/abs(Dc[c[1]]) < TOL*1e2
@test (abs(Dc[c[2]]-λ[2])/abs(Dc[c[2]]) < TOL && abs(Dc[c[3]]-λ[3])/abs(Dc[c[3]]) < TOL*1e2) ||
      (abs(Dc[c[3]]-λ[2])/abs(Dc[c[3]]) < TOL && abs(Dc[c[2]]-λ[3])/abs(Dc[c[2]]) < TOL*1e2)




println("\nTesting in only Complex64 precision (32 bit in real and 32 bit in imaginary)")
nep = nep_gallery("pep0",30)
TOL = 1e-4;
λ,u = jd(Complex64, nep, tol=TOL, maxit=25, Neig = 2, displaylevel = 1, v0=ones(size(nep,1)))
println(" Resnorm of computed solution: ",compute_resnorm(nep,λ[1],u[:,1]))
println(" Smallest eigevalue found: \n λ: ",λ)
Dc,Vc = polyeig(Complex64,nep,DefaultEigSolver)
c = sortperm(abs.(Dc))
println(" 4 smallest eigenvalues according to the absolute values: \n ", Dc[c[1:4]])

@test norm(compute_Mlincomb(nep,λ[1],u[:,1])) < TOL
@test norm(compute_Mlincomb(nep,λ[2],u[:,2])) < TOL
@test  abs(Dc[c[1]]-λ[1])/abs(Dc[c[1]]) < TOL*50
@test (abs(Dc[c[2]]-λ[2])/abs(Dc[c[2]]) < TOL*50) || (abs(Dc[c[3]]-λ[2])/abs(Dc[c[3]]) < TOL*50) #Complex conjugate eigenvalues



println("\nTesting SG as inner solver")
nep = nep_gallery("real_quadratic")
nep = SPMF_NEP(get_Av(nep), get_fv(nep))
TOL = 1e-10;
# Also test that a warning is issued
@test_warn "maxit = 60 is larger than size of NEP = 4. Setting maxit = size(nep,1)" λ,u=jd(Float64, nep, tol=TOL, maxit=60, displaylevel = 1, projtype = :Galerkin, inner_solver_method = NEPSolver.SGIterInnerSolver, v0=ones(size(nep,1)))
λ,u=jd(Float64, nep, tol=TOL, maxit=4, displaylevel = 1, projtype = :Galerkin, inner_solver_method = NEPSolver.SGIterInnerSolver, v0=ones(size(nep,1)))
λ = λ[1]
u = vec(u)
println(" Resnorm of computed solution: ",compute_resnorm(nep,λ,u))
println(" Smallest eigevalue found: \n λ: ",λ)

@test norm(compute_Mlincomb(nep,λ,u)) < TOL



println("\nTesting IAR as projected solver")
nep = nep_gallery("dep0_sparse",40)
TOL = 1e-10;
λ,u = jd(Complex128, nep, tol=TOL, maxit=30, displaylevel = 1, inner_solver_method = NEPSolver.IARInnerSolver, v0=ones(size(nep,1)))
λ = λ[1]
u = vec(u)
println(" Resnorm of computed solution: ",compute_resnorm(nep,λ,u))
println(" Smallest eigevalue found: \n λ: ",λ)

@test norm(compute_Mlincomb(nep,λ,u)) < TOL



println("\nTesting convergence before starting")
λ,u=jd(nep, tol=TOL, maxit=25, Neig=1, displaylevel=1, v0=ones(size(nep,1)), λ=λ, v0=u)
λ = λ[1]
u = vec(u)
@test norm(compute_Mlincomb(nep,λ,u)) < TOL



println("\nTesting errors thrown")
nep = nep_gallery("pep0",4)
# SG requires Galerkin projection type to keep Hermitian
@test_throws ErrorException λ,u=jd(Float64, nep, tol=TOL, maxit=4, projtype = :PetrovGalerkin, inner_solver_method = NEPSolver.SGIterInnerSolver, v0=ones(size(nep,1)))
# An undefined projection type
@test_throws ErrorException λ,u=jd(nep, tol=TOL, maxit=4, projtype = :MYNOTDEFINED, v0=ones(size(nep,1)))
# Too many required eigenvalues, will not converge and hence throw an exception
@test_throws NEPCore.NoConvergenceException λ,u=jd(nep, tol=TOL, maxit=4, Neig=1000, v0=ones(size(nep,1)))


end
