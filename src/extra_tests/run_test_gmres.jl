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


solver = GMRESLinSolver(nep, λ, maxiter=1, restart=200)
println("  type = ", typeof(solver))

b = rand(Complex128, 200)

x, conv_history = lin_solve(solver, b)

A = compute_Mder(nep ,λ, 0)
x2 = A\b

println("GMRES relative residual norm = ", norm(compute_Mlincomb(nep,λ,x, a=[1]) - b)/norm(b))
println("Relative error = ", norm(x-x2)/norm(x2))

