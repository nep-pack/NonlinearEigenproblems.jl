push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using LinSolvers
using Gallery

println("Test Res-Inv")


nep=nep_gallery("dep0")


λ=NaN;
x=NaN
c=ones(size(nep,1));
try
    λ,x =resinv(nep,λ=-0.3,displaylevel=1,linsolvercreator=backslash_linsolvercreator,c=c);
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ=e.λ
    x=e.v
end
println(compute_resnorm(nep,λ,x))

println("λ:",λ);
