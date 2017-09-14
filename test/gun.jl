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

# Test compute_Mlincomb:
λ=150^2+2im;
V=randn(n,2);
z1=compute_Mlincomb(nep1,λ,V, [1.0], 1)
z2=compute_Mlincomb(nep_org,λ,V, [1.0], 1)
@test norm(z1-z2)<sqrt(eps())

# Check that two steps of quasinewton always gives the same result
λ_org=0
try
    quasinewton(nep_org,maxit=2,λ=150^2+1im,v=ones(n),displaylevel=1)
catch e
    λ_org=e.λ
end


λ1=0
try
    quasinewton(nep1,maxit=2,λ=150^2+1im,v=ones(n),displaylevel=1)
catch e
    λ1=e.λ
end

@test abs(λ1-λ_org)<sqrt(eps())


