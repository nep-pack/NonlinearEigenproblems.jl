workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory

using PyPlot
using PyCall

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver
using Gallery



#explicit import needed for overloading functions from packages
import NEPCore.compute_Mlincomb

n=100;
A0=rand(n,n);   A1=rand(n,n);
nep=SPMF_NEP([eye(n), A0, A1],[λ->-λ,λ->eye(λ),λ->expm(-λ)])


compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)
(λ,Q,err)=iar_chebyshev(nep,σ=0,γ=1,Neig=20,v=ones(4),displaylevel=0,maxit=100,tol=eps()*100,check_error_every=1,displaylevel=1)

errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end

m=size(err,1);
for i=1:m
    semilogy(3:3:m, err[3:3:m,i],color="black")
end
ylim(1e-16,1e1)
