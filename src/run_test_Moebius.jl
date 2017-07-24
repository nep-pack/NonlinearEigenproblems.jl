workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using gplot_module


nep=nep_gallery("pep0",100)
#nep=nep_gallery("dep0",100)

#a,b,c,d=rand(4,1)+rand(4,1)*im
#a,b,c,d=[1,0,0,1]
a,b,c,d=[1,-im,1,im]    # Cayley




nep1=MobiusTransformNEP(nep,a=a,b=b,c=c,d=d);

target=rand()+rand()*im
s,x=newton(nep1,λ=target,v=ones(size(nep1,1)),displaylevel=1,armijo_factor=0.5,maxit=30)
println("Resnorm Mobius transformed NEP:",norm(compute_Mlincomb(nep1,s,x)))
println("Mobius transformed NEP eig s=",s);

s=(a*s+b)/(c*s+d)

println("Resnorm original NEP:", norm(compute_Mlincomb(nep,s,x)))
println("Original NEP eig s=",s);


println("Runing IAR");
maxit=100
s,X,err=iar(nep1,displaylevel=1,Neig=20,maxit=maxit,σ=0)
gcloseall()
for i=1:maxit
    gsemilogy(1:1:maxit, err[1:1:maxit,i])
end
show_gplot()

s=(a*s+b)./(c*s+d)

errormeasure=default_errmeasure(nep);
for i=1:length(s)
 println("Original nep. Eigenvalue =",s[i]," residual = ",errormeasure(s[i],X[:,i]))
end
