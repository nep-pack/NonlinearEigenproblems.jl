push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
using GalleryRand
using Test
using LinearAlgebra


@testset "GalleryRand" begin

    @testset "Floats" begin
        @testset "Uniformity" begin
        @info "Testing uniformity of floats"
            my_rng = MSWS_RNG()
            n = 100000
            A = zeros(Float64,n)
            for i = 1:n
                A[i] = gen_rng_float(my_rng)
            end
            sort!(A);

            x = A;
            expr_cdf = cumsum(ones(n))./n

            function uniform_cdf(x)
                if x < 0; return 0.0;
                elseif (x >= 0) && (x <= 1); return (x-0)/(1-0);
                else return 1.0; end;
            end
            true_cdf = zeros(Float64,n)
            for i = 1:n
                true_cdf[i] = uniform_cdf(x[i])
            end

            @test maximum(abs.(expr_cdf-true_cdf)) < sqrt(1/n)
            @test A[1] < 10/n
            @test 1-A[n] < 10/n
        end
    end

    @testset "Matrices" begin
        @testset "Dense" begin
        @info "Testing interfaces for dense matrices"
            my_rng = MSWS_RNG()
            A = gen_rng_mat(my_rng, 40, 50)
            @test size(A) == (40,50)
            @test eltype(A) == Float64
            @test !issparse(A)
        end

        @testset "Sparse" begin
        @info "Testing interfaces for sparse matrices"
            my_rng = MSWS_RNG()
            A = gen_rng_spmat(my_rng, 40, 50, 0.4)
            @test size(A) == (40,50)
            @test eltype(A) == Float64
            @test issparse(A)
        end
    end
end
