# Run tests on nep_gallery (not tested elsewhere)

# Intended to be run from nep-pack/ directory or nep-pack/test directory
push!(LOAD_PATH, string(@__DIR__, "/../src"))

using NEPCore
using NEPTypes
using LinSolvers
using NEPSolver
using Gallery
using LinearAlgebra
using Test


@testset "Newton iterations" begin
    println("Testing sine");
    nep=nep_gallery("sine")
    tol=1e-10;
    λ,v=quasinewton(Float64,nep,λ=-4.2,v=ones(size(nep,1)),tol=tol,
                    displaylevel=1,armijo_factor=0.5,armijo_max=20)
    normalize!(v);
    @test norm(compute_Mlincomb(nep,λ,v))<tol*100
end
