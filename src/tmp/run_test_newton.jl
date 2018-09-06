push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinSolvers

#=println("Running Newton on random dep")
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



nep=nep_gallery("dep0")



println("Running Newton on random dep with Float32 arithmetic only")
λ,x =newton(Float32,nep,displaylevel=1,v=[-1,2,-0.1,0.3,0],λ=1.74,maxit=30);

println("Running Newton on random dep with Float32 arithmetic only (with armijo)")
λ,x =newton(Float32,nep,displaylevel=1,v=[-1,2,-0.1,0.3,0],λ=1.74,maxit=30,armijo_factor=0.9);

println("Solution:",(λ,x))
println("Resnorm:",compute_resnorm(nep,λ,x))

nep=nep_gallery("dep0", 50)

println("Running ResInv ")
λ,x =resinv(nep, displaylevel=1,v=ones(size(nep,1)),λ=-3+0.2im,maxit=200);
println("Running ResInv with Armijo")
λ,x =resinv(nep, displaylevel=1,v=ones(size(nep,1)),λ=-3+0.2im,maxit=200,armijo_factor=0.9);

println("Running ResInv with GMRES as solver")

#SETUP GMRES AS LINEAR SOLVER WITH SOME GIVEN PARAMETERS
linsolverkwargs = ((:maxiter,50), (:restart,50))
function my_gmres_linsolvercreator(nep::NEP, λ)
    return gmres_linsolvercreator(nep, λ, linsolverkwargs)
end

λ,x =resinv(nep, displaylevel=1, linsolvercreator=my_gmres_linsolvercreator);

println("Resnorm:",compute_resnorm(nep,λ,x))



nep=nep_gallery("dep0", 50)


println("Running Aug-Newton")
λ,x =augnewton(nep, displaylevel=1,v=ones(50),λ=-1,maxit=50);
println("Resnorm:",compute_resnorm(nep,λ,x), " eig:",λ)

println("Running Aug-Newton (with armijo)")
λ,x =augnewton(nep, displaylevel=1,v=ones(50),λ=-1,maxit=50,armijo_factor=0.5);
println("Resnorm:",compute_resnorm(nep,λ,x), " eig:",λ)


println("Running quasinewton (without armijo)")
λ,x =quasinewton(Float64,nep, λ=3.5, displaylevel=1,v=ones(size(nep,1),1)[:]);

println("Resnorm:",compute_resnorm(nep,λ,x)/opnorm(x), " eig:",λ)

println("Running quasinewton (with armijo)")
λ,x =quasinewton(Float64,nep, λ=3.5, displaylevel=1,v=ones(size(nep,1),1)[:],armijo_factor=0.9,armijo_max=10);=#


nep = nep_gallery("dep_double");
println("Running implicit det")
n=size(nep,1);
λ,x = implicitdet(nep, λ=9im, v=ones(n), displaylevel=1);
λ,x = newton(nep, λ=9im, v=ones(n), displaylevel=1,maxit=100);
println(λ)
println("Newton QR")
λ,x,y = newtonqr(nep, λ=9im, v=ones(n), displaylevel=1 );

#λ,x = (nep, λ=-0.36, v=ones(n), displaylevel=1 );
#println("Resnorm:",compute_resnorm(nep,λ,x), " eig:",λ)
