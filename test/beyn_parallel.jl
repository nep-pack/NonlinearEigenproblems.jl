# Run parallelization tests on Beyns contour integral method

# Start julia with julia -p XX where XX is the number of cores/processes
# on your computer. On eight.math.kth.se you typically want to do
# julia -p 20 

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
using BenchmarkTools


nep=nep_gallery("dep0",500)

println("Contour integral with quadgk")
tt1=@timev λ1,v1=contour_beyn(nep,displaylevel=0,radius=0.3,k=3,quad_method=:quadgk,tolerance=1e-5,displaylevel=0)
display(tt1)

println("Contour integral with quadg")
tt2=@timev λ2,v2=contour_beyn(nep,displaylevel=0,radius=0.3,k=3,quad_method=:quadg,tolerance=1e-5)
display(tt2)

println("Contour integral with quadg_parallel")

tt3=@timev λ3,v3=contour_beyn(nep,displaylevel=0,radius=0.3,k=3,quad_method=:quadg_parallel,tolerance=1e-5)
display(tt3)


println("Contour integral with quadgk (again for verification)")
tt1b=@timev λ1,v1=contour_beyn(nep,displaylevel=0,radius=0.3,k=3,quad_method=:quadgk,tolerance=1e-5,displaylevel=0)
display(tt1b)


