workspace()
push!(LOAD_PATH, pwd()) # looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using LinSolvers
using Gallery

nep1 = nep_gallery("pep0");
#nep = nep_gallery("dep0_sparse",200);
#
#t=1;
#minusop= S-> -S
#oneop= S -> eye(S)
#expmop= S -> expm(full(-t*S))
#fi=[minusop, oneop, expmop];
#
#nep1=SPMF_NEP([speye(size(nep,1)),nep.A[1],nep.A[2]],fi)
#

print("############### Testing Non-Linear Arnoldi with the PEP: \"pep0\" #################\n");

v=ones(size(nep1,1));
D,X = nlar(nep1,maxit=200,λ=1,nev=10,displaylevel=1,v=v);

#λ,Q,err = iar(pep,maxit=100,Neig=2,σ=0.0,γ=2);
print(norm(compute_Mlincomb(nep1,D[1],X[:,1])))
