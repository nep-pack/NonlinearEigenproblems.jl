workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using PyPlot
using PyCall

n=4;
A0=[0.3000   -0.6000         0    0.4000
   -0.3000    0.4000   -0.8000    1.9000
    0.1000   -1.6000   -1.3000         0
   -1.4000   -0.9000    0.2000    0.9000];

A1=[0.8000    0.2000   -1.3000   -0.3000
   -1.1000    0.9000    1.2000    0.5000
    0.5000    0.2000   -1.6000   -1.3000
    0.7000    0.4000   -0.4000         0];


# define the quadratic DEP
# IDEA: should we define PDEP=polynomial DEP?
nep=SPMF_NEP([eye(4), A0, A1],[λ->-λ^2,λ->eye(λ),λ->expm(-λ)])


λ,Q,err,V,H = iar_chebyshev(nep, maxit=30,Neig=10,σ=0.0,γ=1,displaylevel=1,check_error_every=1);
errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end

m=size(err,1);
for i=1:m
    semilogy(1:1:m, err[1:1:m,i],color="black")
end
