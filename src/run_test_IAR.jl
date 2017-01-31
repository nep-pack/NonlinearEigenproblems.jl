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
nep=nep_gallery("dep0")

function compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))
    return compute_Mlincomb_from_Mder(nep,λ,V,a)
end


λ,Q,err = iar(nep,maxit=100,Neig=1);

for i=1:100
 semilogy(1:100, err[1:100,i], color="red", linewidth=2.0, linestyle="--")
end

errormeasure=default_errmeasure(nep);
for i=1:length(λ)
 println("Residual ",i," eigenpair = ",errormeasure(λ[i],Q[:,i]))
end

