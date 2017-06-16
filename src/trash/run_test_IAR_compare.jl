workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery


# explicit import needed for overloading
# functions from packages
import NEPCore.compute_Mlincomb

println("Load dep0")
nep=nep_gallery("dep0",2000)
#nep=nep_gallery("pep0");

function compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))
    return compute_Mlincomb_from_Mder(nep,λ,V,a)
end

m=3;

tic()
λ,Q,err = iar(nep,maxit=m,Neig=m,σ=2.0,γ=0.5);
toc()

tic()
λ,Q,err = iar_old(nep,maxit=m,Neig=m,σ=2.0,γ=0.5);
toc()
