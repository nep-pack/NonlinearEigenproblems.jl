push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
using GalleryRand
using Test
using LinearAlgebra
using Plots



@testset "GalleryRand" begin

    @testset "Integers" begin
    @info "Testing uniformity of integers"
        reset_rng!()
        n = 1000000;
        upp = 97;
        prob = zeros(Int64,upp); # Avoid zero-indexing
        for i = 1:n
            prob[gen_rng_int(upp)+1] += 1;
        end
        @test iszero(sum(iszero.(prob)))
        true_prob = n/upp;
        @test abs(minimum(prob)-true_prob)/abs(true_prob) < 100*1/sqrt(n)
        @test abs(maximum(prob)-true_prob)/abs(true_prob) < 100*1/sqrt(n)
    end

    @testset "Float interval" begin
    @info "Testing uniformity of floats in an interval"
        reset_rng!()
        n = 1000000;
        A = zeros(Float64,n);
        low = 3.0; upp = 12.0;
        for i = 1:n
            A[i] = gen_rng_float(low, upp);
        end
        sort!(A);

        x = A;
        expr_cdf = cumsum(ones(n))./n;

        function uniform_cdf(x)
            if x < low; return 0.0;
            elseif (x >= low) && (x <= upp); return (x-low)/(upp-low);
            else return 1.0; end;
        end
        true_cdf = zeros(Float64,n)
        for i = 1:n
            true_cdf[i] = uniform_cdf(x[i]);
        end

        @test maximum(abs.(expr_cdf-true_cdf)) < sqrt(1/n)
        @test A[1]-low < 10/n
        @test upp-A[n] < 10/n
    end
end
