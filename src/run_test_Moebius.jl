workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery




nep=nep_gallery("dep0",10)

#a,b,c,d=rand(4,1)+rand(4,1)*im
#a,b,c,d=[1,0,0,1]
a,b,c,d=[1,-im,1,im]    # Cayley
#a,b,c,d=[im,im,-1,1]    # Cayley

#s=(d*s-b)/(-c*s+a);


nep1=MobiusTransformNEP(nep,a=a,b=b,c=c,d=d);

import NEPCore.compute_Mlincomb
compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM(nep,λ,V,a)
compute_Mlincomb(nep::MobiusTransformNEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_from_MM(nep,λ,V,a)

import NEPCore.compute_Mder
compute_Mder(nep::DEP,λ::Number,i::Integer=0)=compute_Mder_from_MM(nep,λ,i)
compute_Mder(nep::MobiusTransformNEP,λ::Number,i::Integer=0)=compute_Mder_from_MM(nep,λ,i)

target=rand()+rand()*im
s,x=newton(nep1,λ=target,v=ones(size(nep1,1)),displaylevel=1,armijo_factor=0.5,maxit=30)
println("Resnorm Mobius transformed NEP:",norm(compute_Mlincomb(nep1,s,x)))
println("Mobius transformed NEP eig s=",s);


#s=(d*s-b)/(-c*s+a);
s=(a*s+b)/(c*s+d);

println("Resnorm original NEP:", norm(compute_Mlincomb(nep,s,x)))
println("Original NEP eig s=",s);




println("Runing IAR");
s,X=iar(nep1,displaylevel=1,Neig=3,σ=0)
s=(a*s+b)./(c*s+d);

for i=1:size(s,1)
    println("Orgnep λ:",s[i]," resnorm orgnep: ",compute_resnorm(nep,s[i],X[:,i]))
end
