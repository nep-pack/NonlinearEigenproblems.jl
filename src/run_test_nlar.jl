workspace()
push!(LOAD_PATH, pwd()) # looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using LinSolvers
using Gallery

pep = nep_gallery("pep0");

print("############### Testing Non-Linear Arnoldi with the PEP: \"pep0\" #################\n");

D,X = nlar(pep,maxit=100,nev=10,displaylevel=1);
