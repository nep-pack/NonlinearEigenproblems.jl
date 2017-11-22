workspace()
push!(LOAD_PATH, string(@__DIR__, "/../src"))

using NEPSolver
using NEPCore
using NEPTypes
using LinSolvers
using Gallery

using Base.Test



#  A Polynomial eigenvalue problem
nep = nep_gallery("pep0",100)

jdtest=@testset "Jacobi–Davidson" begin

λ,u =jd(nep, tol=1e-10, maxit=80, Neig = 5, displaylevel=1)
println("\n Resnorm of computed solution: ",compute_resnorm(nep,λ[1],u[:,1]))
println("\n Smallest eigevalues found: \n λ: ",λ)
Dc,Vc = polyeig(nep,DefaultEigSolver)
c = sortperm(abs.(Dc))
println("\n 5 smallest eigenvalues according to the absolute values: \n", Dc[c[1:5]])

@test norm(compute_Mlincomb(nep,λ[1],u[:,1])) < 1e-10
@test norm(compute_Mlincomb(nep,λ[2],u[:,2])) < 1e-10
@test norm(compute_Mlincomb(nep,λ[3],u[:,3])) < 1e-10
@test norm(compute_Mlincomb(nep,λ[4],u[:,4])) < 1e-10

@test (abs(Dc[c[1]]-λ[1]) < 1e-10 && abs(Dc[c[2]]-λ[2]) < 1e-10) || (abs(Dc[c[2]]-λ[1]) < 1e-10 && abs(Dc[c[1]]-λ[2]) < 1e-10) #Complex conjugate eigenvalues
@test  abs(Dc[c[3]]-λ[3]) < 1e-10
@test (abs(Dc[c[4]]-λ[4]) < 1e-10 && abs(Dc[c[5]]-λ[5]) < 1e-10) || (abs(Dc[c[5]]-λ[4]) < 1e-10 && abs(Dc[c[4]]-λ[5]) < 1e-10) #Complex conjugate eigenvalues



println("\n Testing in only Complex64 precision (32 bit in real and 32 bit in imaginary)")
λ,u =jd(Complex64, nep, tol=1e-4, maxit=80, Neig = 2, displaylevel = 1)
println("\n Resnorm of computed solution: ",compute_resnorm(nep,λ[1],u[:,1]))
println("\n Smallest eigevalue found: \n λ: ",λ)
Dc,Vc = polyeig(Complex64,nep,DefaultEigSolver)
c = sortperm(abs.(Dc))
println("\n 5 smallest eigenvalues according to the absolute values: \n", Dc[c[1:5]])

@test norm(compute_Mlincomb(nep,λ[1],u[:,1])) < 1e-4
@test norm(compute_Mlincomb(nep,λ[2],u[:,2])) < 1e-4
@test (abs(Dc[c[1]]-λ[1]) < 1e-4 && abs(Dc[c[2]]-λ[2]) < 1e-4) || (abs(Dc[c[2]]-λ[1]) < 1e-4 && abs(Dc[c[1]]-λ[2]) < 1e-4) #Complex conjugate eigenvalues



println("\n Testing an (almost) non-PEP")
nep2 = nep_gallery("real_quadratic")
nep2 = SPMF_NEP(get_Av(nep2), get_fv(nep2))
λ,u =jd(Float64, nep2, tol=1e-10, maxit=80, displaylevel = 1)
λ = λ[1]
u = vec(u)
println("\n Resnorm of computed solution: ",compute_resnorm(nep2,λ,u))
println("\n Smallest eigevalue found: \n λ: ",λ)

@test norm(compute_Mlincomb(nep2,λ,u)) < 1e-10



println("\n Testing IAR as projected solver")
# nep_org_3 = nep_gallery(NLEVP_NEP,"gun")
# nep3 = nlevp_make_native(nep_org_3)
# nep3 = nep_gallery("dep0")
nep3 = nep_gallery("pep0",300)
f = function(::Type{T_orig_nep}, T, proj_nep, N) where {T_orig_nep <: ProjectableNEP}
    λ,Q,err,V = iar(T, proj_nep, Neig=1, maxit=100, tol=1e-11)
    return λ[1], Q[:,1]
    # resinv(T, proj_nep, tol = 1e-11, maxit = 100, displaylevel = 1)
end

λ,u =jd(Complex128, nep3, tol=1e-10, maxit=80, displaylevel = 1, proj_eig_solver = f)
λ = λ[1]
u = vec(u)
println("\n Resnorm of computed solution: ",compute_resnorm(nep3,λ,u))
println("\n Smallest eigevalue found: \n λ: ",λ)

@test norm(compute_Mlincomb(nep3,λ,u)) < 1e-10

end
