workspace()
push!(LOAD_PATH, pwd()) # looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using LinSolvers
using Gallery

pep = nep_gallery("pep0");

print("############### Testing Non-Linear Arnoldi with the PEP: \"pep0\" #################\n");

D,X = nlar(pep,maxit=200,λ=1,nev=10,displaylevel=1);

#λ,Q,err = iar(pep,maxit=100,Neig=2,σ=0.0,γ=2);
print(norm(compute_Mlincomb(pep,D[1],X[:,1])))
