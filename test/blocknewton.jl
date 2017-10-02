# Run tests on block Newton

# Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, string(@__DIR__, "/../src"))
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using Base.Test

nep=nep_gallery("dep0",3);

#V=randn(size(nep,1),2);
V=eye(size(nep,1),2);
#V,=qr(V,thin=true)
S0=zeros(2,2);
S,V=blocknewton(nep,S=S0,X=V,displaylevel=1)

λv=eigvals(S);
@test minimum(svdvals(compute_Mder(nep,λv[1])))<sqrt(eps())
@test minimum(svdvals(compute_Mder(nep,λv[2])))<sqrt(eps())
