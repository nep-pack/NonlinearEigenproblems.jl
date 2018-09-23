module TestUtils

using InteractiveUtils
using BenchmarkTools
using Statistics
using Printf

export @bench
export @onlybench
export is_test_script
export set_benchmark_duration_seconds
export enable_benchmark
export displaylevel
export set_displaylevel

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
        :(benchit(() -> $(esc(ex))))
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

end
