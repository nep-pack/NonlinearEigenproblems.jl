#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using Gallery_old

println("Testing Augmented Newton")

A=[ones(12,12),eye(12)]
nep=DEP(A,[0,3]);


#CC=compute_Mder(nep,λ,2)
#V=randn(12,2)
#z=compute_Mlincomb(nep,2+1im,V)
#S=randn(3,3);
#Z=compute_MM(nep,S,V)
## Load a delay eigenvalue problem
#nep=nep_gallery("dep0")
#
#
λ=NaN;
x=NaN
try
    λ,x =aug_newton2(nep,displaylevel=1);
catch e
    # Only catch NoConvergence 
    isa(e, NoConvergenceException) || rethrow(e)  
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ=e.λ
    x=e.v
end
println(compute_resnorm(nep,λ,x))
