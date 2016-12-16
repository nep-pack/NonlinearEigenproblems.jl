#  A Polynomial eigenvalue problem
workspace()
push!(LOAD_PATH, pwd())	# look for modules in the current directory
using NEPSolver
using NEPCore
#using PyPlot

n=200; # mat size
p=10; # Poly degree

A=Array{Float64}(n,n,p)
srand(0)
for i=1:p
   A[:,:,i]=randn(n,n);
end

# Create the derivative function for the PEP
function PEP_Md(λ,i=0)
    # Only workds for i=0 or i=1
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
ev=[Inf]
myerrmeasure=function (λ,v)
    e=nep.relresnorm(λ,v)
    global ev=[ev e]
    return e
end

# 

λ,x =newton_raphson(nep,maxit=30,errmeasure=myerrmeasure);
#λ,x =res_inv(nep,maxit=100,errmeasure=myerrmeasure);	# decrease n, e.g., n=20

println("Resnorm of computed solution: ",norm(nep.Md(λ)*x))

## Plot the iteration history
#Pkg.add("PyPlot")
#using PyPlot
#semilogy(ev[2:end])
#ylabel("resnorm");
#xlabel("iteration");
