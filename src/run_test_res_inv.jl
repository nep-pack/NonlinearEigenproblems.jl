#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore

println("hello world from NEP-pack")
println("Hello again")
n=5;
srand(0) # reset the random seed
A0=randn(n,n);
A1=randn(n,n);
I=eye(n,n);
tau=1;

#nep.n=n;
function DEP_Md(λ,i=0)
    # short-hand definition of functions
    if (i==0)
       return -λ*I+A0+A1*exp(-tau*λ)
    else
       return -I-tau*A1*exp(-tau*λ)
    end
end

nep=NEP(n,DEP_Md);
#println(nep)

λ=NaN;
x=NaN
try
    λ,x =res_inv(nep,displaylevel=1);
catch e
    # Only catch NoConvergence 
    isa(e, NoConvergenceException) || rethrow(e)  
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ=e.λ
    x=e.v
end
println(nep.resnorm(λ,x))







