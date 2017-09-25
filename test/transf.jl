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





# Check that PEP transformations correctly transform coefficients
pep0=nep_gallery("pep0",10);
σ=1+0.3im;
α=3;
pep1=shift_and_scale(pep0,shift=σ,scale=α)
λ,v= quasinewton(pep0,λ=1+1im);
norm(compute_Mlincomb(pep0, λ,v))
λv,V=polyeig(pep1);
@test minimum(abs.(λv-(λ-σ)/α))<eps()*100

# Check that real PEP with real transformation is still real
σ=3;
α=1
pep2=shift_and_scale(pep0,shift=σ,scale=α)
@test !isa(Complex,eltype(pep2.A[1])) # Preserve realness






