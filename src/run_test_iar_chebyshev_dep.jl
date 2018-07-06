# workspace()
# push!(LOAD_PATH, pwd())	# looks for modules in the current directory
#
# #using PyPlot
# #using PyCall
#
# using NEPCore
# using NEPTypes
# using LinSolvers
# using NEPSolver
# using Gallery

import NEPSolver.iar_chebyshev;
include("../src/method_iar_chebyshev.jl");


#explicit import needed for overloading functions from packages
import NEPCore.compute_Mlincomb


nep=nep_gallery("dep0_tridiag",1000)
#nep=nep_gallery("dep0");


#compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)
#(λ,Q,err)=iar(nep,σ=0,γ=1,Neig=6,v=ones(4),displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)
σ=1; γ=0.4; n=size(nep,1)
nep=DEP([broadcast(*,nep.A,exp.(-nep.tauv*σ)/γ); [-σ/γ*eye(n)] ],[nep.tauv*γ; 0])
(λ,Q,err)=iar_chebyshev(nep,σ=-1,γ=2,Neig=6,v=ones(4),displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)
#nep=shift_and_scale(nep,shift=-1,scale=2);
#(λ,Q,err)=iar_chebyshev(nep,σ=0,γ=1,Neig=6,v=ones(4),displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)

errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end
#
# m=size(err,1);
# for i=1:m
#     semilogy(3:3:m, err[3:3:m,i],color="black")
# end
# ylim(1e-16,1e1)
