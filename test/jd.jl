if !isdefined(:global_modules_loaded)
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

@test_throws ErrorException jd()


@testset "Betcke-Voss" begin
println("\n\nTest: Betcke-Voss")


println("\nTesting a PEP")
nep = nep_gallery("pep0",60)
TOL = 1e-11;
λ,u = jd_betcke(nep, tol=TOL, maxit=55, Neig = 3, displaylevel=1, v0=ones(size(nep,1)))
println(" Smallest eigevalues found: \n λ: ",λ)
Dc,Vc = polyeig(nep,DefaultEigSolver)
c = sortperm(abs.(Dc))
println(" 6 smallest eigenvalues according to the absolute values: \n ", Dc[c[1:6]])

# Test residuals
@test norm(compute_Mlincomb(nep,λ[1],u[:,1])) < TOL
@test norm(compute_Mlincomb(nep,λ[2],u[:,2])) < TOL
@test norm(compute_Mlincomb(nep,λ[3],u[:,3])) < TOL



println("\nTesting SG as inner solver")
nep = nep_gallery("real_quadratic")
nep = SPMF_NEP(get_Av(nep), get_fv(nep))
TOL = 1e-10;
# Also test that a warning is issued
λ,u=jd_betcke(Float64, nep, tol=TOL, maxit=4, displaylevel = 1, projtype = :Galerkin, inner_solver_method = NEPSolver.SGIterInnerSolver, v0=ones(size(nep,1)))
λ = λ[1]
u = vec(u)
println(" Resnorm of computed solution: ",compute_resnorm(nep,λ,u))
println(" Smallest eigevalue found: \n λ: ",λ)

@test norm(compute_Mlincomb(nep,λ,u)) < TOL



println("\nTesting IAR as projected solver")
nep = nep_gallery("dep0_sparse",40)
TOL = 1e-10;
λ,u = jd_betcke(Complex128, nep, tol=TOL, maxit=30, displaylevel = 1, inner_solver_method = NEPSolver.IARInnerSolver, v0=ones(size(nep,1)))
λ = λ[1]
u = vec(u)
println(" Resnorm of computed solution: ",compute_resnorm(nep,λ,u))
println(" Smallest eigevalue found: \n λ: ",λ)

@test norm(compute_Mlincomb(nep,λ,u)) < TOL



println("\nTesting convergence before starting")
λ,u=jd_betcke(nep, tol=TOL, maxit=25, Neig=1, displaylevel=1, v0=ones(size(nep,1)), λ=λ, v0=u)
λ = λ[1]
u = vec(u)
@test norm(compute_Mlincomb(nep,λ,u)) < TOL



println("\nTesting errors thrown")
nep = nep_gallery("pep0",4)
# Throw error if iterating more than the size of the NEP
@test_throws ErrorException λ,u=jd_betcke(nep, tol=TOL, maxit=60, displaylevel = 1, v0=ones(size(nep,1)))
# SG requires Galerkin projection type to keep Hermitian
@test_throws ErrorException λ,u=jd_betcke(Float64, nep, tol=TOL, maxit=4, projtype = :PetrovGalerkin, inner_solver_method = NEPSolver.SGIterInnerSolver, v0=ones(size(nep,1)))
# An undefined projection type
@test_throws ErrorException λ,u=jd_betcke(nep, tol=TOL, maxit=4, projtype = :MYNOTDEFINED, v0=ones(size(nep,1)))
# Too many required eigenvalues, will not converge and hence throw an exception
@test_throws NEPCore.NoConvergenceException λ,u=jd_betcke(nep, tol=TOL, maxit=4, Neig=1000, v0=ones(size(nep,1)))

end



@testset "Effenberger" begin
println("\n\nTest: Effenberger")

TOL = 1e-8
nep = nep_gallery("pep0",200)

λ, u = @time jd_effenberger(nep, Neig=3, displaylevel=1, tol=TOL, maxit=100, v0=ones(Complex128,size(nep,1)))
println(" Eigevalues found: \n λ: ",λ)
@test norm(compute_Mlincomb(nep,λ[1],u[:,1])) < TOL
@test norm(compute_Mlincomb(nep,λ[2],u[:,2])) < TOL
@test norm(compute_Mlincomb(nep,λ[3],u[:,3])) < TOL


# λ, u = @time jd_effenberger(nep, Neig=3, displaylevel=1, tol=TOL, maxit=100, inner_solver_method = NEPSolver.ContourBeynInnerSolver, v0=ones(Complex128,size(nep,1)))
# println(" Eigevalues found: \n λ: ",λ)
# @test norm(compute_Mlincomb(nep,λ[1],u[:,1])) < TOL
# @test norm(compute_Mlincomb(nep,λ[2],u[:,2])) < TOL
# @test norm(compute_Mlincomb(nep,λ[3],u[:,3])) < TOL


println("\nTesting errors thrown")
# Throw error if iterating more than the size of the NEP
@test_throws ErrorException λ, u = jd_effenberger(nep, tol=TOL, maxit=(size(nep,1)+1), v0=ones(size(nep,1)))
# SG not possible with Effenberger
@test_throws ErrorException λ, u = jd_effenberger(nep, tol=TOL, maxit=100, inner_solver_method = NEPSolver.SGIterInnerSolver, v0=ones(size(nep,1)))
# Too many required eigenvalues, will not converge and hence throw an exception
@test_throws NEPCore.NoConvergenceException λ, u = jd_effenberger(nep, Neig=1000, tol=TOL, maxit=5, v0=ones(size(nep,1)))


end

end
