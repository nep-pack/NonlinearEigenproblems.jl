include("../src/method_tiar.jl");
nep=nep_gallery("dep_symm_double",1000)
tiar(nep,σ=-1,γ=2;neigs=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
