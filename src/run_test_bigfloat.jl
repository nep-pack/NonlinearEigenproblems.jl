#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using Gallery


println("Test NEP-test with BigFloat")
A0=ones(BigFloat,4,4)-eye(BigFloat,4,4)
u=Array{BigFloat,1}(1:4); v=u-2;
A1=u*v';
A2=eye(BigFloat,4,4);
A2[2,1]=BigFloat(pi);


nep=PEP([A0,A1,A2]);

# start values need to be bigfloats as well
v0=ones(BigFloat,4);
λ0=BigFloat(0)
# 
λ,v=aug_newton(nep,v=v0,λ=λ0,
               tolerance=BigFloat(1e-60),displaylevel=1)

norm(compute_Mlincomb(nep,λ,v))







