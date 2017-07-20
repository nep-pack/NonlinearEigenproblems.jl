workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery




nep=nep_gallery("dep0")

σ=0.8;
γ=0.4;
nep1=ShiftScaleNEP(nep,shift=σ,scale=γ);


s,x=augnewton(nep1,λ=0,v=ones(size(nep1,1)),displaylevel=1)

println("Resnorm shift and scaled NEP:",norm(compute_Mlincomb(nep1,s,x)))
println("Shift and scaled eigval s=",s);

λ=γ*s+σ;

println("Resnorm original NEP:", norm(compute_Mlincomb(nep,λ,x)))

println("Orgnep λ=",λ);
