#  A Polynomial eigenvalue problem
workspace()
push!(LOAD_PATH, pwd())	# look for modules in the current directory
using NEPSolver
using NEPCore
using Gallery


nep=nep_gallery("pep0")

# Saving the errors in an array
ev=zeros(0)
myerrmeasure=function (λ,v)
    e=nep.relresnorm(λ,v)
    global ev=[ev;e]
    return e
end

# 

println("Running Newton Raphson")
λ,x =newton_raphson(nep,maxit=30,errmeasure=myerrmeasure,displaylevel=1);

λ_exact=λ
ev2=zeros(0)
abserrmeasure=function (λ,v)
    e=abs(λ-λ_exact) # Error measure: abs error in λ 
    global ev2=[ev2;nep.resnorm(λ,v)]
    return e
end
println("Running residual inverse iteration")


λ0=round(λ_exact*10)/10; # Start value
x0=round(x*10)/10;       # Start vector
λ,x =res_inv(nep,λ=λ0,v=x0,
             errmeasure=abserrmeasure,displaylevel=1);	

println("Resnorm of computed solution: ",norm(nep.Md(λ)*x))


## Slow first time it
#Pkg.add("Winston")
#using Winston
#semilogy(ev)
#hold(true)
#semilogy(ev2,"r")
#
## Corresponding code in PyPlot
 #Pkg.add("PyPlot")
#using PyPlot
#semilogy(ev[2:end])
#semilogy(ev2[:])
#ylabel("resnorm");
#xlabel("iteration");
