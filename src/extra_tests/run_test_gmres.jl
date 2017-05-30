workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using Suppressor # Pkg.add("Suppressor");
@suppress_err begin
using NEPSolver
using NEPCore
using NEPTypes

using LinSolvers
using Gallery




println("===========================")
println("||  This is GMRES-test   ||")
println("===========================")



nep = nep_gallery("pep0_sparse_003")
λ = rand(Complex128)

import LinSolvers.GMRESLinSolver
gmres_kwargs = ((:maxiter,200), (:restart,200), (:log,true))
function my_gmres_linsolvercreator(nep::NEP, λ)
    return gmres_linsolvercreator(nep, λ, gmres_kwargs)
end

solver = my_gmres_linsolvercreator(nep, λ)
println("  type = ", typeof(solver))

b = rand(Complex128, 200)

x,conv = lin_solve(solver, b)

println(conv[:resnorm])

A = compute_Mder(nep ,λ, 0)
x2 = A\b

println("GMRES relative residual norm = ", norm(compute_Mlincomb(nep, λ, x, a=[1]) - b)/norm(b))
println("Relative error = ", norm(x-x2)/norm(x2))

end
