#using NonlinearEigenproblemsTest
using NonlinearEigenproblems
using NonlinearEigenproblems.RKHelper
using Random
using SparseArrays
using Test

@bench @testset "cached_lin_solver" begin
    dep = nep_gallery("dep0")
    Random.seed!(0)
    y = randn(dep.n)

    cache = nothing

    function solve_and_verify(v, add_to_cache = true)
        x = solve(cache, v, y, add_to_cache)
        @test compute_Mder(dep, v) * x ≈ y
    end

    test_size(cache_size) = @test length(cache.solver) == cache_size

    function run_tests(cache_creator)
        cache = cache_creator(); test_size(0)

        solve_and_verify(2.5, false); test_size(0)
        solve_and_verify(2.5); test_size(1)
        solve_and_verify(2.5); test_size(1)
        solve_and_verify(3.5); test_size(2)

        cache = cache_creator(); test_size(0)

        solve_and_verify(2.5); test_size(1)
    end

    run_tests(() -> LinSolverCache(Float64, dep, DefaultLinSolverCreator()))
    run_tests(() -> LinSolverCache(Float64, dep, BackslashLinSolverCreator()))
end
