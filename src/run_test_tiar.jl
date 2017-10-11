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
# import NEPCore.compute_Mlincomb

println("Load dep0")
nep=nep_gallery("dep0",100)

function compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))
    return compute_Mlincomb_from_MM!(nep,λ,V,a)
end

m=200;   p=1;
λ,Q,err,Z = tiar(nep,maxit=m,Neig=4,σ=2.0,γ=3,check_error_every=p,displaylevel=1);

errormeasure=default_errmeasure(nep);
for i=1:length(λ)
 println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end

m=size(err,1);
for i=1:m
    semilogy(3:3:m, err[3:3:m,i],color="black")
end
ylim(ymax=100)

#λ,Q,err,Z = tiar(nep,maxit=m,Neig=4,σ=λ[1],γ=3,check_error_every=p,displaylevel=1);
