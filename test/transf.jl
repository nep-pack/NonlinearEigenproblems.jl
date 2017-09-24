# Run tests for the NEP transformations

# Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, pwd()*"/src")
push!(LOAD_PATH, pwd()*"/src/gallery_extra")
push!(LOAD_PATH, pwd()*"/src/gallery_extra/waveguide")

using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using Base.Test



orgnep=nep_gallery("dep0");

# Test the shift by solving the orgnep and
# check with residual of transformed nep.
σ=-0.3+0.1im;
nep1=shift_and_scale(orgnep,shift=σ);
orgλ,orgv=augnewton(orgnep)
@test norm(compute_Mlincomb(nep1,orgλ-σ,orgv))<100*eps()
λ1,v1=quasinewton(nep1)
@test abs(λ1+σ-orgλ)<eps()*100 # check that we get the same eigvals




# Test shift and scaling
σ=-0.4+0.01im; α=0.5
nep2=shift_and_scale(orgnep,shift=σ,scale=α);
λ2,v2=quasinewton(nep2)
@test abs((α*λ2+σ)-orgλ)<eps()*100

