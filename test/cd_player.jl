# Run tests on the fiber problem in NLEVP (bessel function nonlinearity)

using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test
using LinearAlgebra

push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
using GalleryNLEVP

@testset "NLEVP cd_player" begin
    nep_org = nep_gallery(NLEVP_NEP,"cd_player")
    nep_nat = nep_gallery("nlevp_native_cd_player")
    n = size(nep_org,1)

    @bench @testset "compare to native" begin
        for seed = [0, 45]
            Random.seed!(seed) # reset the random seed
            λ = rand(ComplexF64)
            for d = [0,1,2,3,4]
                M_org = compute_Mder(nep_org, λ, d)
                M_nat = compute_Mder(nep_nat, λ, d)
                @test opnorm(M_org - M_nat)/(opnorm(M_org)+eps()) < 1e-14
            end
        end

    end

end
