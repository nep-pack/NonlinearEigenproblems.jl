workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory

using PyPlot
using PyCall

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver
using Gallery


#explicit import needed for overloading functions from packages
import NEPCore.compute_Mlincomb
#nep=nep_gallery("dep0_tridiag",1000)


srand(0)
n=100;
d=3;
A = Array{Array{Float64}}(d+1)
for j=0:d
    A[j+1]=rand(n,n)
end
nep=PEP(A)


compute_Mlincomb(nep::PEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)

#(λ,Q,err)=iar(nep,σ=0,γ=1,Neig=10,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)
(λ,Q,err)=iar_chebyshev(nep,σ=0,γ=1,Neig=10,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)


errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end

m=size(err,1);
for i=1:m
    semilogy(3:3:m, err[3:3:m,i],color="black")
end
ylim(1e-16,1e1)
