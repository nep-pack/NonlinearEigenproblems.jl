workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinSolvers


# explicit import needed for overloading
# functions from packages
import NEPCore.compute_Mlincomb

nep=nep_gallery("pep0_sparse_003",500)
#nep=nep_gallery("dep0",100)
#nep=nep_gallery("pep0");

compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)
#
m=50;
λ,Q,err = iar(nep,maxit=m,Neig=5,σ=2.0,γ=3);

errormeasure=default_errmeasure(nep);

# EXTRACT EIGENVECOR
v0=rand(size(nep,1));
println("Residual before extraction = ",errormeasure(λ[1],v0))



v=compute_eigvec_from_eigval_lu(nep,λ[1],default_linsolvercreator);
println("Residual after extraction new = ",errormeasure(λ[1],v))

# example: assume that M0inv already exists, the linsolvercreator is the
# naive function linsolvercreator = (nep, σ) -> M0inv)
M0inv = default_linsolvercreator(nep,λ[1])
v=compute_eigvec_from_eigval_lu(nep,λ[1], (nep, σ) -> M0inv);
println("Residual after extraction = ",errormeasure(λ[1],v))


M0inv = backslash_linsolvercreator(nep,λ[1])
v=compute_eigvec_from_eigval_lu(nep,λ[1], (nep, σ) -> M0inv);
println("Residual after extraction = ",errormeasure(λ[1],v))
