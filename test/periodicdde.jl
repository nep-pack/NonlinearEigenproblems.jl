# Run tests on Periodic DDE

using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test

push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
using GalleryPeriodicDDE

@bench @testset "PeriodicDDE" begin
    nep=nep_gallery(PeriodicDDE_NEP,name="mathieu")
    (λ,v)=quasinewton(ComplexF64,nep,λ=-0.2,v=[1;1])
    # Compare with a reference solution
    @test λ≈(-0.24470143590830754)
end
