




#  A Polynomial eigenvalue problem
workspace()
push!(LOAD_PATH, pwd())	# look for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using LinSolvers
using Gallery



nep = nep_gallery("pep0",100)


λ,u =jd_quad(nep, tol=1e-10, maxit=80);
println("\n Resnorm of computed solution: ",compute_resnorm(nep,λ,u))
println("\n Smallest eigevalue found: \n λ: ",λ)
Dc,Vc = polyeig(nep,DefaultEigSolver);
c = sortperm(abs.(Dc));
println("\n 5 smallest eigenvalues according to the absolute values: \n", Dc[c[1:5]]);



println("\n Testing in only Complex64 precision (32 bit in real and 32 bit in imaginary)")
λ,u =jd_quad(Complex64, nep, tol=1e-4, maxit=80);
println("\n Resnorm of computed solution: ",compute_resnorm(nep,λ,u))
println("\n Smallest eigevalue found: \n λ: ",λ)
Dc,Vc = polyeig(Complex64,nep,DefaultEigSolver);
c = sortperm(abs.(Dc));
println("\n 5 smallest eigenvalues according to the absolute values: \n", Dc[c[1:5]]);



println("\n Testing a non-PEP")
nep2 = nep_gallery("real_quadratic")
λ,u =jd_quad(Complex128, nep2, tol=1e-10, maxit=80);
println("\n Resnorm of computed solution: ",compute_resnorm(nep2,λ,u))
println("\n Smallest eigevalue found: \n λ: ",λ)
