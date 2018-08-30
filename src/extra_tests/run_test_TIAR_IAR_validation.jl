cd("..")
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinearAlgebra

# explicit import needed for overloading
# functions from packages
#import NEPCore.compute_Mlincomb

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
z=rand(100,1);
λ1,Q1,err1 = tiar(nep,maxit=m,Neig=m,σ=2.0,γ=3,v0=z);
λ2,Q2,err2 =  iar(nep,maxit=m,Neig=m,σ=2.0,γ=3,v0=z);

println("Error in the base= ",norm(Q1-Q2))
println("Error in the Ritz values= ",norm(λ1-λ2))
