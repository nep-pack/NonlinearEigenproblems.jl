# Run tests on the fiber problem in NLEVP (bessel function nonlinearity)

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
        for λ = [0.8236475079774124 + 0.9103565379264364im, 0.28604423575660154 + 0.5273506342660073im]
            for d = [0,1,2,3,4]
                M_org = compute_Mder(nep_org, λ, d)
                M_nat = compute_Mder(nep_nat, λ, d)
                @test opnorm(M_org - M_nat)/(opnorm(M_org)+eps()) < 1e-14
            end
        end

    end

end
