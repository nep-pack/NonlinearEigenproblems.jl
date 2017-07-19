workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinSolvers


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


λ,x =newton(Float32,nep,displaylevel=1);

println("Solution:",(λ,x))
println("Resnorm:",compute_resnorm(nep,λ,x))



println("Running ResInv with GMRES as solver")

#SETUP GMRES AS LINEAR SOLVER WITH SOME GIVEN PARAMETERS
linsolverkwargs = ((:maxiter,50), (:restart,50))
function my_gmres_linsolvercreator(nep::NEP, λ)
    return gmres_linsolvercreator(nep, λ, linsolverkwargs)
end

nep=nep_gallery("dep0", 50)
λ,x =resinv(nep, displaylevel=1, linsolvercreator=my_gmres_linsolvercreator);

println("Resnorm:",compute_resnorm(nep,λ,x))



println("Running Aug-Newton")
nep=nep_gallery("dep0", 50)


λ,x =augnewton(nep, displaylevel=1);

println("Resnorm:",compute_resnorm(nep,λ,x), " eig:",λ)


println("Running quasinewton (without armijo)")
λ,x =quasinewton(Float64,nep, λ=3.5, displaylevel=1,v=ones(size(nep,1),1)[:]);

println("Resnorm:",compute_resnorm(nep,λ,x), " eig:",λ)

println("Running quasinewton (with armijo)")
λ,x =quasinewton(Float64,nep, λ=3.5, displaylevel=1,v=ones(size(nep,1),1)[:],armijo_factor=0.9,armijo_max=10);

println("Resnorm:",compute_resnorm(nep,λ,x), " eig:",λ)




