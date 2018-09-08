push!(LOAD_PATH, normpath(@__DIR__, "..")); using TestUtils
using Random
using SparseArrays
using Test

include(normpath(string(@__DIR__), "..", "..", "src", "nleigs", "lusolver.jl"))

@bench @testset "cached_lu_solver" begin
    new_cache() = LUCache(Float64, false)

    cache = new_cache()

    function lusolve_and_verify(fun, v, y)
        x = lusolve(cache, fun, v, y)
        @test fun(v) * x ≈ y
    end

    test_size(cache_size) = @test length(cache.lu) == cache_size

    function run_tests(fun, y)
        cache = new_cache()

        lusolve_and_verify(fun, 2.5, y); test_size(1)
        lusolve_and_verify(fun, 2.5, y); test_size(1)
        lusolve_and_verify(fun, 3.5, y); test_size(2)

        cache = new_cache()

        lusolve_and_verify(fun, 2.5, y); test_size(1)
    end

    Random.seed!(0)
    n = 20
    y = collect(1.0:n)

    test_size(0)

    # sparse matrix
    A = sprandn(n, n, .1)
    fun = v -> A + im * v * I
    run_tests(fun, y)

    # full matrix
    A = randn(n, n)
    fun = v -> A + im * v * I
    run_tests(fun, y)
end
