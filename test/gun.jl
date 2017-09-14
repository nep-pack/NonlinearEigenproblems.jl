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




# Check comput_MM is correct (by comparing against diagonalization of S)
# by using the identity MM(V,S)=MM(V,W*D*inv(W))=MM(V*W,D)*inv(W) and MM(V1,D)
# for diagonal D can be computed with compute_Mlincomb
V=randn(n,3)+randn(n,3);
S=randn(3,3);

# Native compute_MM
N1=compute_MM(nep1,S,V);

# compute same quantity with diagonalization of S and nep_org

# Diagonalize S
(d,W)=eig(S);
D=diagm(d);
V1=V*W;
# 
N2=hcat(compute_Mlincomb(nep_org,d[1],V1[:,1]),
        compute_Mlincomb(nep_org,d[2],V1[:,2]),
        compute_Mlincomb(nep_org,d[3],V1[:,3]))*inv(W)
@test norm(N1-N2)<sqrt(eps())



V=randn(n,2);

# Test that two version of native SPMF nep produce equal results
z1a=compute_Mlincomb_from_Mder(nep1,10+10im,V,[1 1]);
z1b=compute_Mlincomb_from_MM(nep1,10+10im,V,[1 1]);

@test norm(z1a-z1b)<eps()*100

# Test that it is equal to the MATLAB-version
z_org=compute_Mlincomb(nep1,10+10im,V,a=[1 1])

@test norm(z_org-z1a)<sqrt(eps())
@test norm(z_org-z2a)<sqrt(eps())


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


