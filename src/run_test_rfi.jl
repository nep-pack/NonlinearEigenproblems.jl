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

#####################################################################################
println("\nRunning two-sided RFI on pep0")
nep=nep_gallery("pep0");
At = [nep.A[1]',nep.A[2]',nep.A[3]'];
nept=PEP(At);
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
println(compute_resnorm(nep,λ,x))


#####################################################################################
println("\nRunning two-sided RFI on real_quadratic")
nep=nep_gallery("real_quadratic");
At = [nep.A[1]',nep.A[2]',nep.A[3]'];
nept=PEP(At);
try
    λ,x =rfi(nep,nept,displaylevel=1,λ=-50.0+0im);
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

#####################################################################################
println("\nRunning two-sided RFI on dep_double")
nep=nep_gallery("dep_double");
At = [nep.A[1]',nep.A[2]'];
nept=DEP(At,nep.tauv);
try
    λ,x =rfi(nep,nept,displaylevel=1, λ=2.7994400031644808e-8 + 9.424777903529886im);
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