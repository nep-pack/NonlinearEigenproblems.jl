#  Tests for the projected problem
workspace()
push!(LOAD_PATH, pwd())	# look for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
#using Winston # For plotting

nep=nep_gallery("pep0");

## Saving the errors in an array
ev=zeros(0)
myerrmeasure=function (λ,v)
    e=compute_resnorm(nep,λ,v)
    global ev=[ev;e]
    return e
end
#
#

println("Running Newton Raphson")
λ,x =newton(nep,maxit=30,errmeasure=myerrmeasure,
            displaylevel=1);
#
λ_exact=λ
ev2=zeros(0)

pnep=Proj_NEP(nep);
V=randn(size(nep,1),2)
Q,R=qr(hcat(V,x)) # Make the eigenspace a part of the projection subspace
set_projectmatrices!(pnep,Q,Q);
λ1,z1=newton(pnep,λ=(λ_exact+0.001))

x1=Q*z1; x1=x1/x1[1];

# should be small since the eigenvector is in the subspace
println("Difference of solution from projected problem:", norm(x/x[1]-x1))

