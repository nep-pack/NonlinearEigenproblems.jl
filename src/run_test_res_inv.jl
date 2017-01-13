#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore_old
using Gallery_old

println("Test Res-Inv")


nep=nep_gallery("dep0")


λ=NaN;
x=NaN
try
    λ,x =res_inv(nep,displaylevel=1);
catch e
    # Only catch NoConvergence 
    isa(e, NoConvergenceException) || rethrow(e)  
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ=e.λ
    x=e.v
end
println(nep.resnorm(λ,x))







