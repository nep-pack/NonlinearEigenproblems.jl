# Run tests on block Newton

# Intended to be run from nep-pack/ directory or nep-pack/test directory
workspace()
push!(LOAD_PATH, string(@__DIR__, "/../src"))
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using GalleryNLEVP
using Base.Test

#@testset "sgiter" begin

    nep = nep_gallery("real_quadratic")
#    nep = nep_gallery("dep_distributed")
#    nep = nep_gallery("pep0_sym")
#    nep = nep_gallery(NLEVP_NEP,"hadeler")
    j = 1;
    λ,v = sgiter(Float64,
                 nep,
                 j,
                 λ = -100,
                 displaylevel = 1,
                 maxit = 100)
#    @test norm(compute_Mlincomb(nep,λ,v)) < eps(Float64)*100
    println(norm(compute_Mlincomb(nep,λ,v)))
    println(λ)

#end
