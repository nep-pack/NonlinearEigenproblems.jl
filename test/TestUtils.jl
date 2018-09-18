module TestUtils

using InteractiveUtils
using BenchmarkTools
using Statistics
using Printf
using HTTP
using Dates
using TimeZones
using JSON

export @bench
export @onlybench
export save_benchmark
export compare_benchmarks
export is_test_script
export set_benchmark_duration_seconds
export enable_benchmark

benchmark_duration_seconds = 1.0
set_benchmark_duration_seconds(duration) = global benchmark_duration_seconds = duration

run_benchmark = false
enable_benchmark(enabled = true) = global run_benchmark = enabled

const SEPARATOR = "_~/" # this string should not occur in benchmark names

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

function is_test_script(file::AbstractString)
    if occursin(r"(?i)\.jl$", file)
        src = read(file, String)

        pos = 1
        while pos <= length(src)
            expr, pos = Meta.parse(src, pos)
            contains_test_macro(expr) && return true
        end
    end
    return false
end

function contains_test_macro(expr::Expr)
    if expr.head == :macrocall && in(expr.args[1], [Symbol("@test") Symbol("@testset")])
        return true
    end
    return any(e -> contains_test_macro(e), filter(a -> isa(a, Expr), expr.args))
end

mutable struct AggregatedTrial
    name
    key
    trial
end

"Hierarchically aggregate benchmark trials"
function aggregate_trials(aggregated_trials, test, chain)
    chain = [chain; test]
    if haskey(benchmark_results, test)
        trial = minimum(benchmark_results[test])
        for i in 1:length(chain)
            t = chain[i]
            name = " "^((i-1)*2) * t.description
            key = i == 1 ? "" : join([x.description for x in chain[2:i]], SEPARATOR)
            aggr = get!(aggregated_trials, t, AggregatedTrial(name, key, zero(trial)))
            aggr.trial += trial
        end
    end
    foreach(t -> aggregate_trials(aggregated_trials, t, chain), test.results)
    return aggregated_trials
end

struct hierarchical_benchmark
    trial
    indented_name
    children

    hierarchical_benchmark(trial, name) = new(trial, name, Dict())
end

recursive_mapreduce(f, op, rec, root) =
    mapreduce(x -> recursive_mapreduce(f, op, rec, x), op, rec(root); init = f(root))

"Compare two JSON file benchmarks, and print the results to stdout"
function compare_benchmarks(baseline_file, alternative_file)
    base_dict = load_benchmark(baseline_file)
    base = base_dict["benchmark"][1]

    alt_dict = load_benchmark(alternative_file)
    alt = alt_dict["benchmark"][1]

    # reconstruct the testset hierarchy
    components(x) = x[1] == "" ? [] : split(x[1], SEPARATOR)
    alts = sort(collect(alt), lt = (a,b) -> length(components(a)) < length(components(b)))
    root = hierarchical_benchmark(alts[1], "All tests")
    for a in alts[2:end]
        comps = components(a)
        parent = reduce((parent, child) -> parent.children[child], comps[1:end-1]; init = root)
        name = " "^(length(comps)*2) * comps[end]
        parent.children[comps[end]] = hierarchical_benchmark(a, name)
    end

    column_width = recursive_mapreduce(x -> length(x.indented_name), max, parent -> values(parent.children), root)
    column_title = "Test name"
    divider = color(WHITE, "─"^(column_width + 39)) * "\n"
    header = color(WHITE, @sprintf("%s%s         Time               Memory",
        column_title, " "^(column_width - length(column_title))))
    @printf("%s%s\n%s", divider, header, divider)

    # hierarchically sort and print benchmark trials
    function print_hierarchy(test)
        key = test.trial[1]
        trial = test.trial[2]

        if haskey(base, key)
            j = judge(trial, base[key])
            time_diff = 100 * (j.ratio.time - 1)
            memory_diff = 100 * (j.ratio.memory - 1)
        else
            time_diff = NaN
            memory_diff = NaN
        end

        @printf("%s%s    %s    %s\n",
            test.indented_name,
            " "^(column_width - length(test.indented_name)),
            color(time_color(time_diff), prettytime(trial.time) * " (" * prettypercent(time_diff) * ")"),
            color(memory_color(memory_diff), prettymemory(trial.memory) * " (" * prettypercent(memory_diff) * ")"))

        sorted_children = sort(
            [x for x in values(test.children)];
            lt = (a,b) -> a.trial[2].time > b.trial[2].time)
        foreach(print_hierarchy, sorted_children)
    end

    print_hierarchy(root)
    print(divider)

    print_details(title, d) = println("\n",
        color(WHITE, title), "\n",
        color(WHITE, "Run:"), " $(d["time"]); $(d["config"])\n",
        color(WHITE, "CPU:"), " $(d["cpu"])\n",
        color(WHITE, "Git:"), " $(d["git"])")
    print_details("Baseline", base_dict)
    print_details("Alternative", alt_dict)
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

function safe_read_command(cmd)
    try
        chomp(read(cmd, String))
    catch
        ""
    end
end

"Save collected benchmark results to a JSON file on disk"
function save_benchmark(test_results, file_name)
    aggregated_trials = aggregate_trials(Dict(), test_results, [])
    benchmark = BenchmarkGroup()
    foreach(t -> benchmark[t.key] = t.trial, values(aggregated_trials))

    dict = Dict()

    dict["time"] = Dates.now(Dates.UTC)
    dict["config"] = @sprintf("seconds=%.1f", benchmark_duration_seconds)
    dict["cpu"] = @sprintf("%s (%s), %d threads, %.0f GB memory",
        Sys.cpu_info()[1].model,
        Sys.CPU_NAME,
        Sys.CPU_THREADS,
        Sys.total_memory() / (1024.0 * 1024 * 1024))

    git_commit = safe_read_command(`git log -1 --format="%h"`)
    git_branch = safe_read_command(`git symbolic-ref --short HEAD`)
    git_details = safe_read_command(`git log -1 --format="%ai, %s"`)
    dict["git"] = @sprintf("%s (%s), %s", git_commit, git_branch, git_details)

    io = IOBuffer()
    BenchmarkTools.save(io, benchmark)
    dict["benchmark"] = JSONText(String(io.data))

    open(file_name, "w") do f
        JSON.print(f, dict)
    end

    println("Benchmark saved to file $(abspath(file_name))")
end

"Load benchmark results from a JSON file on disk"
function load_benchmark(file_name)
    dict = JSON.parsefile(file_name)

    # convert UTC time to formatted local time
    dt = DateTime(dict["time"])
    zdt = astimezone(ZonedDateTime(dt, tz"UTC"), localzone())
    dict["time"] = Dates.format(zdt, "yyyy-mm-dd HH:MM:SS")

    io = IOBuffer()
    JSON.print(io, dict["benchmark"])
    dict["benchmark"] = BenchmarkTools.load(IOBuffer(io.data))

    return dict
end

const WHITE = 15

color(c, text) = @sprintf("\e[38;5;%dm%s\e[0m", c, text)

const TIME_TOLERANCE_PERCENT = [
    ( -10, 10)
    ( -5, 114)
    (  5, 7)
    (  10, 214)
    (Inf, 9)]

time_color(diff) =
    TIME_TOLERANCE_PERCENT[findfirst(t -> (isnan(diff) ? 0 : diff) < t[1], TIME_TOLERANCE_PERCENT)][2]

const MEMORY_TOLERANCE_PERCENT = [
    ( -5, 10)
    ( -1, 114)
    (  1, 7)
    (  5, 214)
    (Inf, 9)]

memory_color(diff) =
    MEMORY_TOLERANCE_PERCENT[findfirst(t -> (isnan(diff) ? 0 : diff) < t[1], MEMORY_TOLERANCE_PERCENT)][2]

prettypercent(p) = lpad(isnan(p) ? "--- " : @sprintf("%+.0f%%", p), 5, " ")

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
    end
    return lpad(str, 8, " ")
end

end
