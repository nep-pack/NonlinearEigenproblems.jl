# Run tests for the waveguide eigenvalue problem

# Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, pwd()*"/src")	
push!(LOAD_PATH, pwd()*"/src/gallery_extra")
push!(LOAD_PATH, pwd()*"/src/gallery_extra/waveguide")	

using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinSolvers

using Base.Test



nep=nep_gallery("waveguide")
n=size(nep,1);
v0=vec(eye(n,1));
λ0=-3-3im
v0=compute_Mder(nep,λ0)\v0; v0=v0/norm(v0);
v0=compute_Mder(nep,λ0)\v0; v0=v0/norm(v0);
v0=compute_Mder(nep,λ0)\v0; v0=v0/norm(v0);
v0=compute_Mder(nep,λ0)\v0; v0=v0/norm(v0);
v0=compute_Mder(nep,λ0)\v0; v0=v0/norm(v0);

n0=norm(compute_Mlincomb(nep,λ0,v0))
myerrmeasure=(λ,v) -> (norm(compute_Mlincomb(nep,λ,v))/(n0*norm(v)))
try 
    λ,v= quasinewton(Complex128,nep,displaylevel=1,λ=λ0,v=v0,
                     errmeasure=myerrmeasure,maxit=100,tol=1e-2
                     )
catch e
    println(e.λ)
end

    
