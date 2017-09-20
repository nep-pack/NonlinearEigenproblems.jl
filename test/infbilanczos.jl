# Test for infinite Bi-Lanczos

workspace()
push!(LOAD_PATH, pwd()*"/src")	
using NEPSolver
using NEPCore
using NEPTypes
using Gallery

# explicit import needed for overloading
# functions from packages
import NEPCore.compute_Mlincomb

println("Load dep0")

nep=nep_gallery("dep0",100)
#nep=nep_gallery("pep0");
println("Transposing dep0");
# Construct the transposed NEP
nept=DEP([nep.A[1]',nep.A[2]'],nep.tauv)



m=20;
λ,Q,err = infbilanczos(nep,nept,maxit=m,Neig=m,σ=2.0,γ=2);
minimum(svdvals(compute_Mder(nep,λ)))

#
##for i=1:m
## semilogy(1:m, err[1:m,i], color="red", linewidth=2.0, linestyle="--")
##end
##
#errormeasure=default_errmeasure(nep);
#for i=1:length(λ)
# println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
#end
#
