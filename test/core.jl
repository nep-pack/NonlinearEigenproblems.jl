# Tests for core functionality


# Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, string(@__DIR__, "/../src"))


using NEPSolver
using NEPCore
using NEPTypes
using Gallery

using Base.Test

nep=nep_gallery("dep0");
n=size(nep,1);

## Mlincomb_tests
V=randn(n,3);
λ=0.3+1im;
z1=compute_Mlincomb(nep,λ,V)
z2=compute_Mlincomb(nep,λ,V,a=[1 1 1])
@test z1==z2
z3=compute_Mder(nep,λ,0)*V[:,1]+ compute_Mder(nep,λ,1)*V[:,2]+ compute_Mder(nep,λ,2)*V[:,3]
@test z1 ≈ z3

## Test Mlincomb starting counting with kth derivative
z2=compute_Mlincomb(nep,λ,V,a=[0 0 1])
z3=compute_Mder(nep,λ,2)*V[:,3]
@test z2≈z3
z4=compute_Mlincomb(nep,λ,V[:,3],[1], 2);
@test z3≈z4
#
