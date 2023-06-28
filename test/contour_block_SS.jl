# Run tests on the block SS contour integral method

using NonlinearEigenproblems
using Test
using LinearAlgebra


@testset "block SS" begin
    nep=nep_gallery("dep0",3)

    # Circle
    λ,V=contour_block_SS(nep,radius=1.0,N=1000,σ=0.1,k=3,K=3);
    @test norm(compute_Mlincomb(nep,λ[1],V[:,1])) < sqrt(eps())

    # Ellipse
    λ,V=contour_block_SS(nep,radius=[1.0,2.0],N=1000,σ=0.1,k=3,K=3);
    @test norm(compute_Mlincomb(nep,λ[1],V[:,1])) < sqrt(eps())

    # JSIAM Mode
    λ,V=contour_block_SS(nep,radius=1.0,N=1000,σ=0.1,k=3,K=4,Shat_mode=:JSIAM);
    @test norm(compute_Mlincomb(nep,λ[1],V[:,1])) < sqrt(eps())

end
