module TestUtils

using InteractiveUtils
using BenchmarkTools
using Statistics
using Printf
using HTTP
using Dates

export @bench
export report_benchmarks

const time_tolerance_percent = 10
const memory_tolerance_percent = 1
const ci_github_repo = "maxbennedich/julia-ci"
const source_github_repo = "nep-pack/NonlinearEigenproblems.jl"

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

    new_run_context = get_run_context()
    new_benchmark = BenchmarkGroup()
    foreach(t -> new_benchmark[t.key] = t.trial, values(aggregated_trials))

    write_file("run_context.txt", new_run_context)
    BenchmarkTools.save("benchmark.json", new_benchmark)

    old_run_context = read_github_file("run_context.txt")
    json = read_github_file("benchmark.json")
    old_benchmark = BenchmarkTools.load(IOBuffer(json))[1]

    open("benchmark.md", "w") do io
        commit_msg = get(ENV, "TRAVIS_COMMIT_MESSAGE", "Unknown commit")
        commit_id = get(ENV, "TRAVIS_COMMIT", "")
        print(io, "### $commit_msg\n")
        isempty(commit_id) || print(io, "commit [$commit_id](https://github.com/$source_github_repo/commit/$commit_id)\n")

        print(io, "```diff\n")
        column_width = maximum([length(test.name) for test in values(aggregated_trials)])
        column_title = "Test name"
        divider = " " * "─"^(column_width + 39) * "\n"
        @printf(io, "%s %s%s         Time               Memory\n%s", divider, column_title, " "^(column_width - length(column_title)), divider)

        # hierarchically sort and print benchmark trials
        function print_trials(test)
            trial = aggregated_trials[test]

            if haskey(old_benchmark, trial.key) && haskey(new_benchmark, trial.key)
                j = judge(new_benchmark[trial.key], old_benchmark[trial.key])
                time_diff = 100 * (j.ratio.time - 1)
                memory_diff = 100 * (j.ratio.memory - 1)
            else
                time_diff = NaN
                memory_diff = NaN
            end

            diff_prefix =
                time_diff >=  time_tolerance_percent || memory_diff >=  memory_tolerance_percent ? '-' :
                time_diff <= -time_tolerance_percent || memory_diff <= -memory_tolerance_percent ? '+' : ' ';

            @printf(io, "%c%s%s    %s (%s)    %s (%s)\n",
                diff_prefix,
                trial.name,
                " "^(column_width - length(trial.name)),
                prettytime(trial.trial.time),
                prettypercent(time_diff),
                prettymemory(trial.trial.memory),
                prettypercent(memory_diff))

            sorted_children = sort(
                [x for x in test.results if haskey(aggregated_trials, x)];
                lt = (a,b) -> isless(aggregated_trials[b].trial.time, aggregated_trials[a].trial.time))
            foreach(print_trials, sorted_children)
        end

        print_trials(test_results)
        print(io, divider)
        print(io, "```\n")

        print(io, "Previous context:\n")
        print(io, "```\n")
        print(io, old_run_context)
        print(io, "```\n")

        print(io, "Current context:\n")
        print(io, "```\n")
        print(io, new_run_context)
        print(io, "```\n")
    end
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

function read_github_file(file)
    req = HTTP.request("GET", "https://raw.githubusercontent.com/$ci_github_repo/master/$file")
    String(req.body)
end

function write_file(file, str)
    open(file, "w") do f
        write(f, str)
    end
end

function get_run_context()
    io = IOBuffer()
    write(io, "UTC time: $(Dates.format(Dates.now(Dates.UTC), "yyyy-mm-dd HH:MM:SS"))\n")
    write(io, "Commit: $(get(ENV, "TRAVIS_COMMIT", "N/A"))\n")
    write(io, "Message: $(get(ENV, "TRAVIS_COMMIT_MESSAGE", "N/A"))\n")
    write(io, "Travis build: $(get(ENV, "TRAVIS_BUILD_NUMBER", "N/A"))\n")
    write(io, "\n")
    versioninfo(io, verbose = true)
    return String(io.data)
end

function prettypercent(p)
    lpad(isnan(p) ? "--- " : @sprintf("%+.0f%%", p), 5, " ")
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
