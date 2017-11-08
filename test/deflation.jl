# Run tests for the deflation
# Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, string(@__DIR__, "/../src"))
push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide"))

using NEPSolver
using NEPCore
using Gallery
using NEPTypes
using Base.Test

nep=nep_gallery("dep0");

(λ,v)=newton(nep);
n=size(nep,1);
S0=reshape([λ],1,1);
V0=reshape(v,n,1);
dnep=effenberger_deflation(nep,S0,V0)

S=randn(2,2);
V=ones(n+1,2);
compute_MM(dnep,S,V)

(λ2,v2)=augnewton(dnep);
@test minimum(svdvals(compute_Mder(nep,λ2)))<eps()*100

