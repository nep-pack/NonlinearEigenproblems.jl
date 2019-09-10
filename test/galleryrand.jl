push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra"))
using GalleryRand
using Test
using LinearAlgebra


@testset "GalleryRand" begin

    @testset "Integers" begin
        @testset "Interval" begin
        @info "Testing integers drawn in correct interval"
            reset_rng!()
            n = 100000;
            low = -7; upp = 22;
            isininterval_a = true;
            isininterval_b = true;
            for i = 1:n
                a = gen_rng_int(low, upp)
                b = gen_rng(low, upp)
                isininterval_a = isininterval_a && (a>=low) && (a<=upp)
                isininterval_b = isininterval_b && (b>=low) && (b<=upp)
            end
            @test isininterval_a
            @test isininterval_b
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
            n = 100000;
            low = -7.1; upp = 22;
            isininterval_a = true;
            isininterval_b = true;
            for i = 1:n
                a = gen_rng_float(low, upp);
                b = gen_rng(low, upp)
                isininterval_a = isininterval_a && (a>=low) && (a<=upp)
                isininterval_b = isininterval_b && (b>=low) && (b<=upp)
            end
            @test isininterval_a
            @test isininterval_b
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

    @testset "Complex" begin
        @testset "Interval" begin
        @info "Testing complex drawn in correct interval"
            reset_rng!()
            n = 100000;
            real_low = -7.1;
            real_upp = 2.9;
            imag_low = -4.0;
            imag_upp = 1;
            a = real_low + imag_upp*1im;
            b = real_upp + imag_low*1im;
            isininterval_c = true;
            isininterval_d = true;
            for i = 1:n
                c = gen_rng_complex(a, b);
                d = gen_rng(a, b);
                isininterval_c = isininterval_c && (real(c)>=real_low) && (real(c)<=real_upp) &&
                                                   (imag(c)>=imag_low) && (imag(c)<=imag_upp)
                isininterval_d = isininterval_d && (real(d)>=real_low) && (real(d)<=real_upp) &&
                                                   (imag(d)>=imag_low) && (imag(d)<=imag_upp)
            end
            @test isininterval_c
            @test isininterval_d
        end
    end

    @testset "Matrices" begin
        @testset "Dense" begin
        @info "Testing interfaces for dense matrices"
            reset_rng!()

            A = gen_rng_mat(4, 5, -3, 0);
            @test size(A) == (4,5)
            @test eltype(A) == Int64
            A = gen_rng_mat(4, -1, 3);
            @test size(A) == (4,4)
            @test eltype(A) == Int64

            A = gen_rng_mat(4, 5, -3.0, 0);
            @test size(A) == (4,5)
            @test eltype(A) == Float64
            A = gen_rng_mat(4, -1, 3.1);
            @test size(A) == (4,4)
            @test eltype(A) == Float64

            A = gen_rng_mat(3, 7, -3+1im, 1);
            @test size(A) == (3,7)
            @test eltype(A) == ComplexF64
            A = gen_rng_mat(3, -3+1im, 1.0-1im);
            @test size(A) == (3,3)
            @test eltype(A) == ComplexF64
        end

        @testset "Sparse" begin
        @info "Testing interfaces for sparse matrices"
            reset_rng!()

            A = gen_rng_spmat(4, 5, 0.4, -3, 0);
            @test size(A) == (4,5)
            @test eltype(A) == Int64
            @test issparse(A)
            A = gen_rng_spmat(4, 0.4, -1, 3);
            @test size(A) == (4,4)
            @test eltype(A) == Int64
            @test issparse(A)

            A = gen_rng_spmat(4, 5, 0.4, -3.0, 0);
            @test size(A) == (4,5)
            @test eltype(A) == Float64
            @test issparse(A)
            A = gen_rng_spmat(4, 0.4, -1, 3.1);
            @test size(A) == (4,4)
            @test eltype(A) == Float64
            @test issparse(A)

            A = gen_rng_spmat(3, 7, 0.4, -3+1im, 1);
            @test size(A) == (3,7)
            @test eltype(A) == ComplexF64
            @test issparse(A)
            A = gen_rng_spmat(3, 0.4, -3+1im, 1.0-1im);
            @test size(A) == (3,3)
            @test eltype(A) == ComplexF64
            @test issparse(A)
        end
    end
end
