workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
push!(LOAD_PATH, ".." )	# looks for modules in the current directory
using NEPSolver
using NEPCore
using Gallery


println("Load dep0")
nep=nep_gallery("dep0")


n=nep.n;

# FOR SMALL k IT WORKS
k=10;
W=rand(n,k);
compute_Mlincomb(nep,0.0,W);

# THE FOLLOWING GIVES OVERFLOW (LARGE k)
k=40;
W=rand(n,k);

compute_Mlincomb(nep,0.0,W);
