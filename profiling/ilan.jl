#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

#import NEPSolver.inner_solve;
include("../src/inner_solver.jl");
#import NEPSolver.ilan;
include("../src/method_ilan.jl");
nep=nep_gallery("dep_symm_double",1000)
ilan(nep,σ=0,γ=1;Neig=5,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)
