#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using Gallery
println("Test MSLP")

#
println("Running MSLP sparse dep (with default eigsolver)")
nep=nep_gallery("dep0_sparse")
λ,x =mslp(nep,displaylevel=1);
println(compute_resnorm(nep,λ,x))

# Buggy julia eigs generates error 
# λ,x =mslp(nep,displaylevel=1,eigsolver="eigs");


println("Running MSLP full dep (with eigsolver=eig)")
nep=nep_gallery("dep0")
λ,x =mslp(nep,displaylevel=1);
println(compute_resnorm(nep,λ,x))












