module TestUtils

using BenchmarkTools
using Statistics
using Printf
using HTTP

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
    warmup(@benchmarkable $test[1] = $f())
    benchmark_results[test[1]] = run((@benchmarkable $f()), seconds=1)
end

mutable struct AggregatedTrial
    name
    key
    trial
end

"Print benchmark results collected by `@bench`"
function report_benchmarks(test_results)
    run_benchmark || return

    aggregated_trials = Dict()

    # hierarchically aggregate benchmark trials
    function aggregate_trials(test, chain)
        chain = [chain; test]
        if haskey(benchmark_results, test)
            trial = minimum(benchmark_results[test])
            for i in 1:length(chain)
                t = chain[i]
                name = " "^((i-1)*2) * t.description
                key = i == 1 ? "" : join([x.description for x in chain[2:i]], " / ")
                aggr = get!(aggregated_trials, t, AggregatedTrial(name, key, zero(trial)))
                aggr.trial += trial
            end
        end
        foreach(t -> aggregate_trials(t, chain), test.results)
    end

    aggregate_trials(test_results, [])

    new_benchmark = BenchmarkGroup()
    foreach(t -> new_benchmark[t.key] = t.trial, values(aggregated_trials))
    BenchmarkTools.save("benchmark.json", new_benchmark)

    json_request = HTTP.request("GET", "https://raw.githubusercontent.com/maxbennedich/julia-ci/master/benchmark.json")
    json = String(json_request.body)
    old_benchmark = BenchmarkTools.load(IOBuffer(json))[1]

    foreach(j -> println("$(j[1]): $(j[2])"), judge(new_benchmark, old_benchmark))
    println("-------------")
    for key in intersect(keys(old_benchmark), keys(new_benchmark))
        j = judge(new_benchmark[key], old_benchmark[key]; time_tolerance = 0.1, memory_tolerance = 0.05)
        if j.time != :invariant || j.memory != :invariant
            println("$key -- $(prettytime(old_benchmark[key].time)) -> $(prettytime(new_benchmark[key].time)) -- $j")
        end
    end

    column_width = maximum([length(test.name) for test in values(aggregated_trials)])
    column_title = "Test name"
    divider = "─"^(column_width + 23) * "\n"
    @printf("%s%s%s     Time       Memory\n%s", divider, column_title, " "^(column_width - length(column_title)), divider)

    # hierarchically sort and print benchmark trials
    function print_trials(test)
        trial = aggregated_trials[test]
        @printf("%s%s    %s    %s\n", trial.name, " "^(column_width - length(trial.name)), prettytime(trial.trial.time), prettymemory(trial.trial.memory))

        sorted_children = sort(
            [x for x in test.results if haskey(aggregated_trials, x)];
            lt = (a,b) -> isless(aggregated_trials[b].trial.time, aggregated_trials[a].trial.time))
        foreach(print_trials, sorted_children)
    end

    print_trials(test_results)
    print(divider)
end

import Base.+
import Base.zero

+(t1::BenchmarkTools.TrialEstimate, t2::BenchmarkTools.TrialEstimate) =
    BenchmarkTools.TrialEstimate(
        t1.params,
        t1.time + t2.time,
        t1.gctime + t2.gctime,
        t1.memory + t2.memory,
        t1.allocs + t2.allocs)

zero(t::BenchmarkTools.TrialEstimate) = BenchmarkTools.TrialEstimate(t.params, 0, 0, 0, 0)

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
