# Run tests on Beyns contour integral method

push!(LOAD_PATH, @__DIR__); using TestUtils
using NonlinearEigenproblems.NEPCore
using NonlinearEigenproblems.NEPSolver
using NonlinearEigenproblems.Gallery
using Test
using LinearAlgebra

@testset "Beyn contour" begin
    nep=nep_gallery("dep0")
    @bench @testset "disk at origin" begin

        λ,v=contour_beyn(nep,displaylevel=1,radius=1,k=2,quad_method=:ptrapz,compute_eigenvectors=true)

        for i = 1:2
            println(i, ":  ", λ[i])
            M=compute_Mder(nep,λ[i])
            minimum(svdvals(M))
            @test minimum(svdvals(M)) < eps()*1000
            @test norm(compute_Mlincomb(nep,λ[i],v[:,i]))/norm(v[:,i]) < eps()*100
        end


        λ,v=contour_beyn(nep,displaylevel=1,radius=1,k=2,quad_method=:ptrapz,N=100)
        M=compute_Mder(nep,λ[1])
        minimum(svdvals(M))
        @test minimum(svdvals(M))<eps()*1000
        @test all(isnan.(v))

    end
    @bench @testset "shifted disk" begin

        λ,v=contour_beyn(nep,displaylevel=1,σ=-0.2,radius=1.5,k=4,quad_method=:ptrapz,compute_eigenvectors=true)

        for i = 1:4
            println(i, ":  ", λ[i])
            M=compute_Mder(nep,λ[i])
            minimum(svdvals(M))
            @test minimum(svdvals(M)) < eps()*1000
            @test norm(compute_Mlincomb(nep,λ[i],v[:,i]))/norm(v[:,i]) < eps()*100
        end

    end

end
