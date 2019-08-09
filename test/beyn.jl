# Run tests on Beyns contour integral method

using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using Test
using LinearAlgebra


@testset "Beyn contour" begin
    nep=nep_gallery("dep0")
    @bench @testset "disk at origin" begin

        λ,v=contour_beyn(nep,logger=displaylevel,radius=1,neigs=2,sanity_check=false)

        for i = 1:2
            @info "$i: $(λ[i])"
            M=compute_Mder(nep,λ[i])
            minimum(svdvals(M))
            @test minimum(svdvals(M)) < eps()*1000
            @test norm(compute_Mlincomb(nep,λ[i],v[:,i]))/norm(v[:,i]) < eps()*500
        end


        λ,v=contour_beyn(nep,logger=displaylevel,radius=1,k=2,N=100,sanity_check=false)
        M=compute_Mder(nep,λ[1])
        minimum(svdvals(M))
        @test minimum(svdvals(M))<eps()*1000

    end
    @bench @testset "shifted disk" begin

        λ,v=contour_beyn(nep,logger=displaylevel,σ=-0.2,radius=1.5,
                         neigs=4,sanity_check=false)

        @test size(λ,1)==4
        for i = 1:4
            @info "$i: $(λ[i])"
            M=compute_Mder(nep,λ[i])
            minimum(svdvals(M))
            @test minimum(svdvals(M)) < eps()*10000
            @test norm(compute_Mlincomb(nep,λ[i],v[:,i]))/norm(v[:,i]) < eps()*10000
        end

    end

end
