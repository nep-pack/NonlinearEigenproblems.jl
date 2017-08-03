workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery




# functions from packages
nep=nep_gallery("dep0_sparse",10);nept=DEP([nep.A[1]',nep.A[2]'],nep.tauv);

#nep=nep_gallery("pep0",50);
#nept=PEP([nep.A[1]',nep.A[2]',nep.A[3]'])

λ,Q,err = iar(nep,maxit=100,Neig=2,σ=1.0,γ=1,displaylevel=0,check_error_every=1);

w=compute_eigvec_from_eigval_lopcg(nep,nept,λ[1]);
errormeasure=default_errmeasure(nep);
println("\n Eigenvalue=",λ[1],"\nResidual after extraction = ",errormeasure(λ[1],w))
