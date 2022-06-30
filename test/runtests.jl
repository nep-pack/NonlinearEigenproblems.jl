using BenchmarkTools, Test, Printf, NonlinearEigenproblems;
################################################################################
# Runs all .jl files in the directory of this script containing a @test or
# @testset macro, except those specified to be excluded in run_all_tests.jl
################################################################################

# This section activates the test project. It's needed for "Pkg.test()" to work for the root project.
push!(LOAD_PATH, "@stdlib") # needed for "using Pkg"
using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate();



displaylevel = 1
set_displaylevel(level) = global displaylevel = level

benchmark_duration_seconds = 1.0
set_benchmark_duration_seconds(duration) = global benchmark_duration_seconds = duration

run_benchmark = false
enable_benchmark(enabled = true) = global run_benchmark = enabled

benchmark_results = Dict()

"""
This macro can be used in front of `@testset` macros. It will only take effect if the
`run_benchmark` global variable is set to `true`. If so, it will do a simple
benchmark of the code within the `@testset`, and store the results within this module,
which can later be reported by calling `save_benchmark`.

Since all benchmark results are aggregated hierarchically, it is advised not to nest
`@bench` expressions, since that will end up double counting elapsed time and memory.
"""
macro bench(ex)
    if run_benchmark
        test = ex.args[end]
        if test.head == :for
            # rewrite "@bench @testset for ... end" as "for ... @bench @testset begin ... end; end"
            rewritten_expr =
                Expr(test.head, test.args[1],
                    Expr(:macrocall, Symbol("@bench"), LineNumberNode(0),
                        Expr(ex.head, ex.args[1:end-1]..., test.args[2])))
            return :($(esc(rewritten_expr)))
        else
            :(benchit(() -> $(esc(ex))))
        end
    else
        :($(esc(ex)))
    end
end

"Like `@bench`, but if `run_benchmark` is false, the `@testset` won't be executed at all."
macro onlybench(ex)
    if run_benchmark
        :(@bench($(esc(ex))))
    end
end

function benchit(f)
    test = Vector(undef, 1)

    warmup(@benchmarkable $test[1] = $f())

    # "Infinite" number of samples, to enforce the configured benchmark duration.
    # There will be a single eval per sample, which at least should be fine for benchmark
    # expressions longer than ~1 μs.
    benchmark_results[test[1]] = run((@benchmarkable $f()), samples = typemax(Int), seconds = benchmark_duration_seconds)
end

"""
`@test` that there are `expected_nr_λ` eigenvalues in the `λ` vector, and that all
`default_errmeasure` norms of the `nep` with corresponding eigenvectors in the `V`
matrix are below the given tolerance.
"""
function verify_lambdas(expected_nr_λ, nep::NEP, λ, V, tol = 1e-5)
    @test length(λ) == expected_nr_λ
    @info "Found $(length(λ)) lambdas:"
    errmeasure=ResidualErrmeasure(nep);
    for i in eachindex(λ)
        nrm = estimate_error(errmeasure,λ[i], V[:, i])
        @test nrm < tol
        @info "λ[$i] = $(λ[i]) (norm = $(@sprintf("%.4g", nrm)))"
    end
end


include("run_all_tests.jl")

test_name_regex = length(ARGS) < 1 ? "" : ARGS[1]
run_all_tests(test_name_regex)
