push!(LOAD_PATH, pwd()) # looks for modules in the current directory
using NonlinearEigenproblems: NEPSolver, NEPCore, NEPTypes, LinSolvers, Gallery
using LinearAlgebra
#using gplot_module


#nep1 = nep_gallery("pep0");
nep = nep_gallery("dep0_sparse",200);

t=1;
minusop= S-> -S
oneop= S -> eye(S)
expmop= S -> exp(full(-t*S))
fi=[minusop, oneop, expmop];

nep1=SPMF_NEP([speye(size(nep,1)),nep.A[1],nep.A[2]],fi)


print("############### Testing Non-Linear Arnoldi with the PEP: \"pep0\" #################\n");

v=ones(size(nep1,1));
#D,X,err_hyst = nlar(nep1,maxit=200,λ=1,nev=10,displaylevel=1,v=v);
maxit=200; nev=4;
D,X,err_hyst = nlar(nep1;maxit=maxit,λ=1,nev=nev,displaylevel=1,v=v);

#λ,Q,err = iar(pep,maxit=100,Neig=2,σ=0.0,γ=2);
print(norm(compute_Mlincomb(nep1,D[1],X[:,1])))

#=
gcloseall()
for i=1:nev-1
  gsemilogy(1:maxit, err_hyst[:,i])
end
show_gplot()=#
