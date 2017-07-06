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
nep=nep_gallery("dep0",100)
#nep=nep_gallery("pep0");
#
function compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))
    #return compute_Mlincomb_from_Mder(nep,λ,V,a)
    return compute_Mlincomb_from_MM!(nep,λ,V,a)
end
#
m=50;
λ,Q,err = tiar(nep,maxit=m,Neig=m,σ=0.0,γ=1);

for i=1:m
 semilogy(1:m, err[1:m,i], color="red", linewidth=2.0, linestyle="--")
end

errormeasure=default_errmeasure(nep);
for i=1:length(λ)
 println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end
