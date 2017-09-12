# Run tests on Newton methods

# Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, pwd()*"/src")	
push!(LOAD_PATH, pwd()*"/../src")	
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinSolvers

using Base.Test

nep=nep_gallery("dep0")

# newton and augnewton are equivalent, therefore I expect them
# to generate identical results
λ1,x1 =newton(nep,displaylevel=0,v=ones(size(nep,1)),λ=0,tolerance=eps()*10);
λ2,x2 =augnewton(nep,displaylevel=0,v=ones(size(nep,1)),λ=0,tolerance=eps()*10);


@test λ1 ≈ λ2
@test x1 ≈ x2
r1=compute_resnorm(nep,λ1,x1)
r2=compute_resnorm(nep,λ2,x2)

@test r1 < eps()*100
@test r2 < eps()*100





