#  A Polynomial eigenvalue problem
workspace()
push!(LOAD_PATH, pwd()) # look for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using LinSolvers
using Gallery

include("method_sg_old.jl")


nep=nep_gallery("real_quadratic");


# Four smallest eigenvalues of the problem "real_quadratic":
# 1.0e+03 *
# -2.051741417993845
# -0.182101627437811
# -0.039344930222838
# -0.004039879577113


λ_init = -100;
λ,x = sg_iteration(nep,tolerance_outer=1e-4,tolerance_inner=1e-4,λ=λ_init,displaylevel=1,eigsolvertype=DefaultEigSolver)
println("Resnorm of computed solution: ",compute_resnorm(nep,λ,x))



λ,x = sg_iteration(Float64, nep,tolerance_outer=1e-8,λ=λ_init,displaylevel=1,eigsolvertype=DefaultEigSolver,maxit=30)
println("Resnorm of computed solution: ",compute_resnorm(nep,λ,x))
