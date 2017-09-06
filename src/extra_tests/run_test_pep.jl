#  A Polynomial eigenvalue problem
workspace()
push!(LOAD_PATH, pwd())	# look for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
#using Winston # For plotting

nep=nep_gallery("pep0");

## Saving the errors in an array
ev=zeros(0)
myerrmeasure=function (λ,v)
    e=compute_resnorm(nep,λ,v)
    global ev=[ev;e]
    return e
end
#
#

println("Running Newton Raphson")
λ,x =newton(nep,maxit=30,errmeasure=myerrmeasure,
            displaylevel=1);
#
λ_exact=λ
ev2=zeros(0)

abserrmeasure=function (λ,v)
    e=abs(λ-λ_exact) # Error measure: abs error in λ 
    global ev2=[ev2;compute_resnorm(nep,λ,v)] # store residual norm
    return e
end
println("Running residual inverse iteration")
#
λ0=round.(λ_exact*10)/10; # Start value
x0=round.(x*10)/10;       # Start vector
λ,x =resinv(nep,λ=λ0,v=x0,displaylevel=1,errmeasure=abserrmeasure);	

println("Resnorm of computed solution: ",compute_resnorm(nep,λ,x))

#
#
## Slow first time it
#Pkg.add("Winston")
#
#semilogy(ev)
#hold(true)
#semilogy(ev2,"r")
#
## Corresponding code in PyPlot
##Pkg.add("PyPlot")
#using PyPlot
#semilogy(ev[2:end])
#semilogy(ev2[:])
#ylabel("resnorm");
#xlabel("iteration");


## Pkg.add("Plots")
using Plots
plotly()
plot(log10(ev[2:end]))
