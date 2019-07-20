# Run tests on the block SS contour integral method

using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test
using LinearAlgebra


@testset "block SS" begin
    nep=nep_gallery("dep0",3)
    λ,V=contour_block_SS(nep,radius=1.0,N=1000,σ=0.1,k=3,K=3);
    @test norm(compute_Mlincomb(nep,λ[1],V[:,1])) < sqrt(eps())
end
