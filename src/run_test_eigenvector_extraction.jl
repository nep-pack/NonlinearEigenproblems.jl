workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using gplot_module

# explicit import needed for overloading
# functions from packages
import NEPCore.compute_Mlincomb

println("Load dep0")
nep=nep_gallery("dep0",100)
#nep=nep_gallery("pep0");

function compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))
    #return compute_Mlincomb_from_Mder(nep,λ,V,a)
    return compute_Mlincomb_from_MM!(nep,λ,V,a)
end
#
m=50;
λ,Q,err = iar(nep,maxit=m,Neig=m,σ=2.0,γ=3);

gcloseall()
mm=50;
for i=1:mm
  gsemilogy(1:m, err[1:m,i])
end
show_gplot()

errormeasure=default_errmeasure(nep);
for i=1:length(λ)
 println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end


# EXTRACT EIGENVECOR
s=λ[1];
v0=rand(size(nep,1));
println("Residual before extraction = ",errormeasure(λ[1],v0))

v=compute_eigvec_from_eigval_old(nep,s;v=v0,tol=sqrt(eps()));
println("Residual after extraction = ",errormeasure(λ[1],v))

v=compute_eigvec_from_eigval(nep,s);
println("Residual after extraction new = ",errormeasure(λ[1],v))
