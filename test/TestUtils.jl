module TestUtils

using BenchmarkTools
using Printf

export @bench
export report_benchmarks

# the TEST_SUITE environment variable can be set by the CI tool
run_benchmark = get(ENV, "TEST_SUITE", "") == "benchmark"

benchmark_results = Dict()

"""
This macro can be used in front of `@testset` macros. It will only take effect if the
`TEST_SUITE` environment variable is set to `benchmark`. If so, it will do a simple
benchmark of the code within the `@testset`, and store the results within this module,
which can later be reported by calling `report_benchmarks`.

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

function benchit(f)
    # Benchmark for a single sample/eval. This will end up running the code twice; once
    # for warmup / JIT compilation, and once for benchmarking.
    test = Vector(undef, 1)
    benchmark_results[test[1]] = @benchmark $test[1] = $f() samples=1 evals=1
end

mutable struct AggregatedResult
    name
    time
    memory
end

"Print benchmark results collected by `@bench`"
function report_benchmarks(test_results)
    run_benchmark || return

    aggregated_results = Dict()

    # hierarchically aggregate benchmark results
    function aggregate_results(test, chain)
        chain = [chain; test]
        if haskey(benchmark_results, test)
            result = benchmark_results[test]
            for i in eachindex(chain)
                t = chain[i]
                aggr = get!(aggregated_results, t, AggregatedResult(" "^((i-1)*2) * t.description, 0.0, 0))
                aggr.time += result.times[1]
                aggr.memory += result.memory
            end
        end
        foreach(t -> aggregate_results(t, chain), test.results)
    end

    aggregate_results(test_results, [])

    column_width = maximum([length(test.name) for test in values(aggregated_results)])
    column_title = "Test name"
    divider = "─"^(column_width + 23) * "\n"
    @printf("%s%s%s     Time       Memory\n%s", divider, column_title, " "^(column_width - length(column_title)), divider)

    # hierarchically sort and print benchmark results
    function print_results(test)
        result = aggregated_results[test]
        @printf("%s%s    %s    %s\n", result.name, " "^(column_width - length(result.name)), prettytime(result.time), prettymemory(result.memory))

        sorted_children = sort(
            [x for x in test.results if haskey(aggregated_results, x)];
            lt = (a,b) -> isless(aggregated_results[b].time, aggregated_results[a].time))
        foreach(print_results, sorted_children)
    end

    print_results(test_results)
    print(divider)
end

# Copied from TimerOutputs.jl
function prettytime(t)
    if t < 1e3
        value, units = t, "ns"
    elseif t < 1e6
        value, units = t / 1e3, "μs"
    elseif t < 1e9
        value, units = t / 1e6, "ms"
    else
        value, units = t / 1e9, "s"
    end

    if round(value) >= 100
        str = string(@sprintf("%.0f ", value), units)
    elseif round(value * 10) >= 100
        str = string(@sprintf("%.1f ", value), units)
    else
        str = string(@sprintf("%.2f ", value), units)
    end
    return lpad(str, 7, " ")
end

# Copied from TimerOutputs.jl
function prettymemory(b)
    if b < 1000
        value = -1
        str = string(round(Int, b), "B")
    elseif b < 1000^2
        value, units = b / 1024, "KiB"
    elseif b < 1000^3
        value, units = b / 1024^2, "MiB"
    else
        value, units = b / 1024^3, "GiB"
    end

    if round(value) >= 100
        str = string(@sprintf("%.0f ", value), units)
    elseif round(value * 10) >= 100
        str = string(@sprintf("%.1f ", value), units)
    elseif value >= 0
        str = string(@sprintf("%.2f ", value), units)
    else
        str = "-"
    end
    return lpad(str, 8, " ")
end

end
