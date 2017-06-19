workspace()
push!(LOAD_PATH, pwd()) # looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinSolvers


println("Running two-sided RFI on random dep")
nep=nep_gallery("dep0")
nept=DEP([nep.A[1]',nep.A[2]'],nep.tauv)

λ=NaN;
x=NaN
try
    λ,x =rfi(nep,nept,displaylevel=1);
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ=e.λ
    x=e.v
end
println(λ)
println("Resnorm:",compute_resnorm(nep,λ,x))

#=println("\nRunning two-sided RFI on dep with double eigenvalue")
nep=nep_gallery("dep_distributed")
try
    λ,x =rfi(nep,nept,displaylevel=1,λ=-0.4+0.9im);
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ=e.λ
    x=e.v
end

println(λ)
println(compute_resnorm(nep,λ,x))=#