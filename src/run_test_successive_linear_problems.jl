#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using Gallery
println("Test MSLP")

nep=nep_gallery("dep0_sparse")

λ=NaN;
x=NaN
try
    λ,x =successive_linear_problems(nep,displaylevel=1,eigsolver="matlab_eigs");
catch e
    # Only catch NoConvergence 
    isa(e, NoConvergenceException) || rethrow(e)  
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ=e.λ
    x=e.v
end
println(nep.resnorm(λ,x))







