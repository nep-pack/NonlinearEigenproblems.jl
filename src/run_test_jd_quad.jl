




#  A Polynomial eigenvalue problem
workspace()
push!(LOAD_PATH, pwd())	# look for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using LinSolvers
using Gallery



n=100; # mat size
	   
A0 = randn(n,n);
A1 = randn(n,n);
A2 = randn(n,n);
A = [A0,A1,A2];
nep=PEP(A);



λ,u =jd_quad(nep,tolerance=1e-10,maxit=80);

println("\n Resnorm of computed solution: ",compute_resnorm(nep,λ,u))
println("\n Smallest eigevalue found: \n λ: ",λ)



Dc,Vc = polyeig(nep,DefaultEigSolver);
c = sortperm(abs(Dc));
println("\n 5 smallest eigenvalues according to the absolute values: \n", Dc[c[1:5]]);


