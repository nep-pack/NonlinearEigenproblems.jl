push!(LOAD_PATH, pwd())	# looks for modules in the current directory
#using Suppressor # Pkg.add("Suppressor");
#@suppress_err begin
using NEPSolver
using NEPCore
using NEPTypes
using LinSolvers
using Gallery
using LinearAlgebra


println("===========================")
println("||  This is GMRES-test   ||")
println("===========================")


nep = nep_gallery("pep0_sparse_003")
λ = rand(ComplexF64)

import LinSolvers.GMRESLinSolver
gmres_kwargs = ((:maxiter,200), (:restart,200), (:log,true))
function my_first_gmres_linsolvercreator(nep::NEP, λ)
    return gmres_linsolvercreator(nep, λ, gmres_kwargs)
end

first_solver = my_first_gmres_linsolvercreator(nep, λ)
println("type = ", typeof(first_solver), "\n")

b = rand(ComplexF64, 200)

x = lin_solve(first_solver, b)


A = compute_Mder(nep ,λ, 0)
x2 = A\b

println("GMRES relative residual norm = ", norm(compute_Mlincomb(nep, λ, x, a=[1]) - b)/norm(b))
println("Relative error = ", norm(x-x2)/norm(x2))



# DO IT AGAIN BUT WITHOUT LOG
gmres_kwargs = ((:maxiter,200), (:restart,200))
function my_second_gmres_linsolvercreator(nep::NEP, λ)
    return gmres_linsolvercreator(nep, λ, gmres_kwargs)
end

second_solver = my_second_gmres_linsolvercreator(nep, λ)

x = lin_solve(second_solver, b)
println("")
println("GMRES relative residual norm = ", norm(compute_Mlincomb(nep, λ, x, a=[1]) - b)/norm(b))
println("Relative error = ", norm(x-x2)/norm(x2))


# DO IT AGAIN BUT WITH VERBOSE AND A PRECONDITIONER
P = lufact(diagm(diag(nep.A[1])+diag(nep.A[2])+diag(nep.A[3])+ones(size(nep,1))))
gmres_kwargs = ((:maxiter,200), (:restart,200), (:log,true), (:verbose,true), (:Pl,P))
function my_third_gmres_linsolvercreator(nep::NEP, λ)
    return gmres_linsolvercreator(nep, λ, gmres_kwargs)
end

third_solver = my_third_gmres_linsolvercreator(nep, λ)

x = lin_solve(third_solver, b)
println("")
println("GMRES relative residual norm = ", norm(compute_Mlincomb(nep, λ, x, a=[1]) - b)/norm(b))
println("Relative error = ", norm(x-x2)/norm(x2))


#end
