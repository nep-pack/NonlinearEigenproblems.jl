#  A Polynomial eigenvalue problem
workspace()
push!(LOAD_PATH, pwd())	# look for modules in the current directory
using NEPSolver
using NEPCore

n=200; # mat size
p=10; # Poly degree

A=Array{Float64}(n,n,p)
srand(0)
for i=1:p
   A[:,:,i]=randn(n,n);
end

# Create the derivative function for the PEP
function PEP_Md(λ,i=0)
    if (i==0)
        M=zeros(n,n)
        for i=1:p
            M+=λ^(i-1)*A[:,:,i]
        end
        return M
    else
        Mp=zeros(n,n)
        for i=2:p
            Mp += (i-1)*(λ^(i-2)*A[:,:,i])
        end
        return Mp
    end
end

nep=NEP(n,PEP_Md);

# Saving the errors in an array
ev=[inf]
myerrmeasure=function (λ,v)
    e=nep.relresnorm(λ,v)
    global ev=[ev e]
    println(e)
    return e
end

# 

λ,x =newtonraphson(nep,maxit=30,errmeasure=myerrmeasure);


println("Resnorm of computed solution: ",norm(nep.Md(λ)*x))

# Plot the iteration history
Pkg.add("PyPlot")
using PyPlot
semilogy(ev[2:end])
ylabel("resnorm")
xlabel("iteration")
