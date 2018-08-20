if !isdefined(:global_modules_loaded)
    using Base.Test
end

include(normpath(string(@__DIR__), "..", "..", "src", "nleigs", "lusolver.jl"))

new_cache() = LUCache(Float64, false)

cache = new_cache()

function lusolve_and_verify(fun, v, y)
    x = lusolve(cache, fun, v, y)
    @test fun(v) * x ≈ y
end

macro testcache(ex, cache_size)
    :(
        $(esc(ex));
        @test length(cache.lu) == $cache_size
    )
end

function run_tests(fun, y)
    global cache = new_cache()

    @testcache lusolve_and_verify(fun, 2.5, y) 1
    @testcache lusolve_and_verify(fun, 2.5, y) 1
    @testcache lusolve_and_verify(fun, 3.5, y) 2

    global cache = new_cache()

    @testcache lusolve_and_verify(fun, 2.5, y) 1
end

srand(0)
n = 20
y = collect(1.0:n)

@testset "cached_lu_solver" begin
    @testcache () 0

    # sparse matrix
    A = sprandn(n, n, .1)
    fun = v -> A + im * v * speye(n, n)
    run_tests(fun, y)

    # full matrix
    A = randn(n, n)
    fun = v -> A + im * v * eye(n, n)
    run_tests(fun, y)
end
