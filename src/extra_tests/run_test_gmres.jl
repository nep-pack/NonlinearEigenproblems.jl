workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
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

gmres_kwargs = ((:maxiter,200), (:restart,200), (:log,true))
solver = GMRESLinSolver(nep, λ; gmres_kwargs...)
println("  type = ", typeof(solver))

b = rand(Complex128, 200)

x,conv = lin_solve(solver, b)

println(conv[:resnorm])

A = compute_Mder(nep ,λ, 0)
x2 = A\b

println("GMRES relative residual norm = ", norm(compute_Mlincomb(nep,λ,x, a=[1]) - b)/norm(b))
println("Relative error = ", norm(x-x2)/norm(x2))

