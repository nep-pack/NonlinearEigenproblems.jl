# Run tests on Beyns contour integral method

using NonlinearEigenproblems
using Test
using Random
using LinearAlgebra


@testset "Beyn contour" begin

    nep=nep_gallery("dep0")
    @bench @testset "disk at origin" begin
        Random.seed!(0);

        λ,v=contour_beyn(nep,logger=displaylevel,radius=1,neigs=1,sanity_check=false)

        for i = 1:size(λ,1)
            @info "$i: $(λ[i])"
            M=compute_Mder(nep,λ[i])
            minimum(svdvals(M))
            @test minimum(svdvals(M)) < eps()*1000
            @test norm(compute_Mlincomb(nep,λ[i],v[:,i]))/norm(v[:,i]) < eps()*500
        end


        λ,v=contour_beyn(nep,logger=displaylevel,radius=0.4,k=2,N=100,sanity_check=false)
        M=compute_Mder(nep,λ[1])
        minimum(svdvals(M))
        @test minimum(svdvals(M))<eps()*1000

    end
    @bench @testset "shifted disk" begin

        λ,v=contour_beyn(nep,logger=displaylevel,σ=0.2,radius=1.0,
                         neigs=4,sanity_check=false)

        @test size(λ,1)==3
        for i = 1:3
            @info "$i: $(λ[i])"
            M=compute_Mder(nep,λ[i])
            minimum(svdvals(M))
            @test minimum(svdvals(M)) < eps()*10000
            @test norm(compute_Mlincomb(nep,λ[i],v[:,i]))/norm(v[:,i]) < eps()*10000
        end

    end

end
