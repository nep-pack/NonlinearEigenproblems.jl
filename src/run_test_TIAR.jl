workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using gplot_module

# explicit import needed for overloading
# functions from packages
#import NEPCore.compute_Mlincomb

println("Load dep0")
nep=nep_gallery("dep0",100)

function compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))
    return compute_Mlincomb_from_MM!(nep,λ,V,a)
end

m=50;   p=3;
λ,Q,err = tiar(nep,maxit=m,Neig=m,σ=2.0,γ=3,p=p);

errormeasure=default_errmeasure(nep);
for i=1:length(λ)
 println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end

gcloseall()
for i=1:m
  gsemilogy(p:p:m, err[p:p:m,i])
end
show_gplot()
