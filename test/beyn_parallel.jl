# Run parallelization tests on Beyns contour integral method

# Start julia with julia -p XX where XX is the number of cores/processes
# on your computer. On eight.math.kth.se you typically want to do
# julia -p 20

# Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, string(@__DIR__, "/../src"))
push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide"))

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver
using Gallery
using IterativeSolvers
using Base.Test
using BenchmarkTools


nep=nep_gallery("dep0",500)

tol=1e-5
println("Contour integral with quadgk")
tt1=@timev λ1,v1=contour_beyn(nep,displaylevel=0,radius=0.3,k=3,quad_method=:quadgk,tol=tol,displaylevel=0)
display(tt1)

println("Contour integral with quadg")
tt2=@timev λ2,v2=contour_beyn(nep,displaylevel=0,radius=0.3,k=3,quad_method=:quadg,tol=tol)
display(tt2)

println("Contour integral with quadg_parallel")

tt3=@timev λ3,v3=contour_beyn(nep,displaylevel=0,radius=0.3,k=3,quad_method=:quadg_parallel,tol=tol)
display(tt3)


println("Contour integral with quadgk (again for verification)")
tt1b=@timev λ1,v1=contour_beyn(nep,displaylevel=0,radius=0.3,k=3,quad_method=:quadgk,tol=1e-5,displaylevel=0)
display(tt1b)

# Check distance between eigenvalue approximations taking
# into account that they are not necessarily ordered
# What we would really want is the Hausdorf distance (but it is
# is not available in the julia core packages)

@test norm(minimum(abs.(λ1*ones(1,3)-ones(3,1)*λ2'),1))<10*tol

@test norm(minimum(abs.(λ3*ones(1,3)-ones(3,1)*λ2'),1))<10*tol
