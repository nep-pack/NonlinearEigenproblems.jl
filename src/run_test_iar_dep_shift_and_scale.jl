workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory

#using PyPlot
#using PyCall

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver
using Gallery

#import NEPSolver.iar_chebyshev;
#include("../src/method_iar_chebyshev.jl");



srand(0) # reset the random seed

n=100;
A1=rand(n,n);
A2=rand(n,n);
A3=rand(n,n);

tau1=0;
tau2=.3;
tau3=2.1;


nep=DEP([A1,A2,A3],[tau1,tau2,tau3])
nep_old=nep;


#explicit import needed for overloading functions from packages
#import NEPCore.compute_Mlincomb
#compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM!(nep,λ,V,a)
σ=1.0; γ=2.0;
(λ,Q,err)=iar(nep,σ=σ,γ=γ,Neig=6,v=ones(4),displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)
# redefine the dep shifted and scaled
mm=length(nep.tauv)
tauv2=zeros(mm+1)
A=Array{Array{Float64,2},1}(mm+1)
for j=1:length(nep.tauv)
	tauv2[j]=nep.tauv[j]*γ;
	A[j]=nep.A[j]*exp(-nep.tauv[j]*σ)/γ
end
A[mm+1]=-σ/γ*eye(n); tauv2[mm+1]=0
nep2=DEP(A,tauv2)

#nep3=DEP([broadcast(*,nep.A,exp.(-nep.tauv*σ)/γ); [-σ/γ*nep.A[1]^0] ],[nep.tauv*γ; 0])
nep3=shift_and_scale(nep,shift=σ,scale=γ);
(μ,Q,err)=iar(nep3,σ=0,γ=1,Neig=6,v=ones(4),displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)
λ2=σ+γ*μ

errormeasure=default_errmeasure(nep);
for i=1:length(λ2)
    println("Eigenvalue=",λ2[i]," residual = ",errormeasure(λ2[i],Q[:,i]))
end
#(λ,Q,err)=iar(nep,σ=0,γ=1,Neig=6,v=ones(4),displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)

#
# errormeasure=default_errmeasure(nep);
# for i=1:length(λ)
#     println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
# end
#
# m=size(err,1);
# for i=1:m
#     semilogy(3:3:m, err[3:3:m,i],color="black")
# end
# ylim(1e-16,1e1)
