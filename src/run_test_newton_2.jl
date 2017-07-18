workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinSolvers


println("Running Newton on random dep")
nep=nep_gallery("dep0")
try
    λ,x =newton(nep,displaylevel=1,maxit=2);
catch e
    println(typeof(e))
end
