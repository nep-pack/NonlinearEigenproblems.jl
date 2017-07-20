workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery




nep=nep_gallery("dep0")

σ=1.5im;
γ=2;
nep1=ShiftScaleNEP(nep,shift=σ,scale=γ);


s,x=newton(nep1,λ=0,v=ones(size(nep1,1)),displaylevel=1,armijo_factor=0.5,maxit=30)

println("Resnorm shift and scaled NEP:",norm(compute_Mlincomb(nep1,s,x)))
println("Shift and scaled eigval s=",s);

λ=γ*s+σ;

println("Resnorm original NEP:", norm(compute_Mlincomb(nep,λ,x)))

println("Orgnep λ=",λ);


println("Runing IAR");
s,X=iar(nep1,displaylevel=1,Neig=3)
λ=γ*s+σ;

for i=1:size(s,1)
    println("Orgnep λ:",λ[i]," resnorm orgnep: ",compute_resnorm(nep,λ[i],X[:,i]))
end
