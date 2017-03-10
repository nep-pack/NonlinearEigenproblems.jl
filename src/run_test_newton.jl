#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery


println("Running Newton on random dep")
nep=nep_gallery("dep0")

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
println("Resnorm:",compute_resnorm(nep,λ,x))

println("\nRunning Newton on dep with double eigenvalue")
nep=nep_gallery("dep_double")
try
    λ,x =newton(nep,displaylevel=1, λ=8.5im, maxit=100);
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



println("Running Newton on random dep with Float32 arithmetic only")
nep=nep_gallery("dep0")


λ,x =newton(Float32,nep,displaylevel=1, tolerance=eps(Float32)*100);

println("Solution:",(λ,x))
println("Resnorm:",compute_resnorm(nep,λ,x))





