workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using PyPlot
using PyCall

n=10
E=rand(2*n)+rand(2*n)*im

nep=nep_gallery("qep_fixed_eig",n,E)

a,b,c,d=[1,-im,1,im]    # Cayley
nep=MobiusTransformNEP(nep,a=a,b=b,c=c,d=d);
#
#f=z->(a*z+b)./(c*z+d)

σ=0; m=100; p=1;
λ,Q,err = iar(nep,maxit=m,Neig=50,σ=σ,γ=3,displaylevel=1,check_error_every=3);
errormeasure=default_errmeasure(nep);
for i=1:length(λ)
    println("Eigenvalue=",λ[i]," residual = ",errormeasure(λ[i],Q[:,i]))
end


for i=1:m
    semilogy(p:p:m, err[p:p:m,i], color="red", linewidth=2.0, linestyle="--")
end
