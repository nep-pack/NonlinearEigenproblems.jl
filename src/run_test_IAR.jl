#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using PyPlot

# explicit import needed for overloading
# functions from packages
import NEPCore.compute_Mlincomb

println("Load dep0")
#nep=nep_gallery("dep0")
nep=nep_gallery("pep0");

function compute_Mlincomb(nep::PEP,λ::Number,V;a=ones(size(V,2)))
    return compute_Mlincomb_from_Mder(nep,λ,V,a)
end

m=30;
λ,Q,err = iar(nep,maxit=m,Neig=m,σ=0.0);

for i=1:m
 semilogy(1:m, err[1:m,i], color="red", linewidth=2.0, linestyle="--")
end

errormeasure=default_errmeasure(nep);
for i=1:length(λ)
 println("Residual ",i," eigenpair = ",errormeasure(λ[i],Q[:,i]))
end

