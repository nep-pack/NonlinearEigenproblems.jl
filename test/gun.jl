# Run tests on Beyns contour integral method 

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



nep_org=nlevp_gallery_import("gun","../nlevp3/");
nep1=nlevp_make_native(nep_org);

n=size(nep1,1);
tol=1e-11;
λ1,v1=quasinewton(nep1,λ=150^2+1im,v=ones(n),displaylevel=1,tolerance=tol,maxit=500);

v1=v1/norm(v1);

@test norm(compute_Mlincomb(nep1,λ1,v1))<tol*100
@test norm(compute_Mder(nep1,λ1)*v1)<tol*100

@test norm(compute_Mlincomb(nep_org,λ1,v1))<tol*100
@test norm(compute_Mder(nep_org,λ1)*v1)<tol*100


