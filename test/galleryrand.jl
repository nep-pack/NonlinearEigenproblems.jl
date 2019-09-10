push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
using GalleryRand
using Test
using LinearAlgebra


@testset "GalleryRand" begin

    @testset "Integers" begin
        @testset "Interval" begin
        @info "Testing integers drawn in correct interval"
            reset_rng!()
            n = 10000;
            low = -7; upp = 22;
            isininterval = true;
            for i = 1:n
                a = gen_rng_int(low, upp)
                isininterval = isininterval && (a>=low) && (a<=upp)
            end
            @test isininterval
        end

        @testset "Uniformity" begin
        @info "Testing uniformity of integers"
            reset_rng!()
            n = 1000000;
            low = 3; upp = 99;
            prob = zeros(Int64, upp); # Avoid zero-indexing
            for i = 1:n
                prob[gen_rng_int(low, upp)] += 1;
            end
            @test sum(iszero.(prob)) == low-1
            prob_int = prob[low:upp];
            true_prob = n/(upp-low+1);
            @test abs(minimum(prob_int)-true_prob)/abs(true_prob) < 100*1/sqrt(n)
            @test abs(maximum(prob_int)-true_prob)/abs(true_prob) < 100*1/sqrt(n)
        end
    end


    @testset "Floats" begin
        @testset "Interval" begin
        @info "Testing floats drawn in correct interval"
            reset_rng!()
            n = 10000;
            low = -7.1; upp = 22;
            isininterval = true;
            for i = 1:n
                a = gen_rng_float(low, upp);
                isininterval = isininterval && (a>=low) && (a<=upp)
            end
            @test isininterval
        end

        @testset "Uniformity" begin
        @info "Testing uniformity of floats"
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
end
