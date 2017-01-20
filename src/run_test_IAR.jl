#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using Gallery


println("Load dep0")
nep=nep_gallery("dep0")


################# NEWTON #########################3

println("Running Newton on random dep")

λ=NaN;
x=NaN
try
    λ,x =newton(nep,displaylevel=1);
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ=e.λ
    x=e.v
end
println(λ)
println(compute_resnorm(nep,λ,x))




D,V = iar(nep);



err=abs(λ-D);
minimum(err);


