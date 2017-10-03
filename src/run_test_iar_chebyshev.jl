workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using PyPlot
using PyCall

# explicit import needed for overloading
# functions from packages
import NEPCore.compute_Mlincomb

nep=nep_gallery("dep0",10)
#nep=nep_gallery("pep0");


compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)

λ,Q,err = iar_chebyshev(nep,maxit=100,Neig=10,σ=2.0,γ=3,displaylevel=1,check_error_every=3);
errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end

m=size(err,1);
for i=1:m
    semilogy(3:3:m, err[3:3:m,i],color="black")
end
