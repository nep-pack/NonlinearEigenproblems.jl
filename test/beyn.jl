# Run tests on Beyns contour integral method 

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

println("Contour integral with σ=0,radius=1,k=2")
λ,v=contour_beyn(nep,displaylevel=1,radius=1,k=2)

println(λ[1])
M=compute_Mder(nep,λ[1])
minimum(svdvals(M))
@test minimum(svdvals(M))<eps()*1000

println(λ[2])
M=compute_Mder(nep,λ[2])
@test minimum(svdvals(M))<eps()*1000


println("Contour integral with σ=-0.2,radius=1.5,k=4")

λ,v=contour_beyn(nep,displaylevel=1,σ=-0.2,radius=1.5,k=4)

println(λ[1])
M=compute_Mder(nep,λ[1])
minimum(svdvals(M))
@test minimum(svdvals(M))<sqrt(eps())

println(λ[2])
M=compute_Mder(nep,λ[2])
@test minimum(svdvals(M))<sqrt(eps())

println(λ[3])
M=compute_Mder(nep,λ[3])
@test minimum(svdvals(M))<sqrt(eps())

println(λ[4])
M=compute_Mder(nep,λ[4])
@test minimum(svdvals(M))<sqrt(eps())
