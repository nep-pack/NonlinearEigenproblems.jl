#using NonlinearEigenproblemsTest
using BenchmarkTools
using Statistics
using Printf
using HTTP
using Dates
using TimeZones
using JSON

include("pretty_formatting.jl")

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

const SEPARATOR = "_~/" # this string should not occur in benchmark names

recursive_mapreduce(f, op, rec, root) =
    mapreduce(x -> recursive_mapreduce(f, op, rec, x), op, rec(root); init = f(root))

mutable struct AggregatedTrial
    name
    key
    trial
end

"Hierarchically aggregate benchmark trials"
function aggregate_trials(aggregated_trials, test, chain)
    chain = [chain; test]
    if haskey(NonlinearEigenproblemsTest.benchmark_results, test)
        trial = minimum(NonlinearEigenproblemsTest.benchmark_results[test])
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

"Save collected benchmark results to a JSON file on disk"
function save_benchmark(test_results, runtime_seconds, file_name)
    aggregated_trials = aggregate_trials(Dict(), test_results, [])
    benchmark = BenchmarkGroup()
    foreach(t -> benchmark[t.key] = t.trial, values(aggregated_trials))

    dict = Dict()

    dict["time"] = Dates.now(Dates.UTC)
    dict["runtime"] = runtime_seconds
    dict["config"] = @sprintf("seconds=%.1f", NonlinearEigenproblemsTest.benchmark_duration_seconds)
    dict["julia-version"] = string(VERSION)
    dict["julia-commit"] = Base.GIT_VERSION_INFO.commit_short
    dict["cpu"] = @sprintf("%s (%s), %d threads, %.0f GB memory",
        Sys.cpu_info()[1].model,
        Sys.CPU_NAME,
        Sys.CPU_THREADS,
        Sys.total_memory() / 1024^3)

    git_commit = safe_read_command(`git log -1 --format="%h"`)
    git_branch = safe_read_command(`git rev-parse --abbrev-ref HEAD`)
    git_branch == "HEAD" && (git_branch = safe_read_command(pipeline(pipeline(`git describe --all --exact-match`, stderr=devnull), `sed 's=.*/=='`)))
    isempty(git_branch) && (git_branch = "N/A")
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
load_benchmark_from_file(file_name) = load_benchmark(JSON.parsefile(file_name))

"Load benchmark results from a JSON string"
load_benchmark_from_string(str) = load_benchmark(JSON.parse(str))

function load_benchmark(dict)
    # convert UTC time to formatted local time
    dt = DateTime(dict["time"])
    zdt = astimezone(ZonedDateTime(dt, tz"UTC"), localzone())
    dict["time"] = Dates.format(zdt, "yyyy-mm-dd HH:MM:SS")

    io = IOBuffer()
    JSON.print(io, dict["benchmark"])
    dict["benchmark"] = BenchmarkTools.load(IOBuffer(io.data))

    return dict
end

function safe_read_command(cmd)
    try
        chomp(read(cmd, String))
    catch
        ""
    end
end

struct hierarchical_benchmark
    trial
    indented_name
    children

    hierarchical_benchmark(trial, name) = new(trial, name, Dict())
end

components(x) = x[1] == "" ? [] : split(x[1], SEPARATOR)

"Reconstruct the testset hierarchy from deserialized benchmark data"
function get_testset_hierarchy(benchmark)
    tests = sort(collect(benchmark), by = x -> length(components(x)))
    root = hierarchical_benchmark(tests[1], "All tests")
    for t in tests[2:end]
        comps = components(t)
        parent = reduce((parent, child) -> parent.children[child], comps[1:end-1]; init = root)
        name = " "^(length(comps)*2) * comps[end]
        parent.children[comps[end]] = hierarchical_benchmark(t, name)
    end
    return root
end

"Hierarchically sort and iterate over benchmark trials"
function iterate_hierarchy(f, test)
    f(test)
    for child in sort([x for x in values(test.children)]; by = t -> -t.trial[2].time)
        iterate_hierarchy(f, child)
    end
end

const WHITE = 15

color(c, text) = @sprintf("\e[38;5;%dm%s\e[0m", c, text)

const TIME_TOLERANCE_PERCENT = [
    (-10, 10)
    ( -5, 114)
    (  5, 7)
    ( 10, 214)
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

function print_table_header(test_hierarchy, time_memory_width, time_padding, memory_padding)
    column_width = recursive_mapreduce(x -> length(x.indented_name), max, parent -> values(parent.children), test_hierarchy)
    column_title = "Test name"
    divider = color(WHITE, "─"^(column_width + time_memory_width)) * "\n"
    header = color(WHITE, @sprintf("%s%sTime%sMemory", column_title,
        " "^(column_width - length(column_title) + time_padding), " "^memory_padding))
    @printf("%s%s\n%s", divider, header, divider)
    return column_width, divider
end

function print_benchmark_details(title, d)
    runstr = d["time"]
    haskey(d, "runtime") && (runstr *= @sprintf("; runtime=%.0fs", d["runtime"]))
    haskey(d, "config") && (runstr *= "; $(d["config"])")
    if haskey(d, "julia-version")
        runstr *= "; Julia $(d["julia-version"])"
        haskey(d, "julia-commit") && (runstr *= " ($(d["julia-commit"]))")
    end

    println("\n",
        color(WHITE, title), "\n",
        color(WHITE, "Run:"), " $runstr\n",
        color(WHITE, "CPU:"), " $(d["cpu"])\n",
        color(WHITE, "Git:"), " $(d["git"])")
end

"Prints a single JSON file benchmark to stdout."
function print_benchmark(baseline_file)
    base_dict = load_benchmark_from_file(baseline_file)
    base = base_dict["benchmark"][1]

    root = get_testset_hierarchy(base)
    column_width, divider = print_table_header(root, 23, 6, 7)

    iterate_hierarchy(root) do test
        trial = test.trial[2]
        @printf("%s%s    %s     %s\n",
            test.indented_name,
            " "^(column_width - length(test.indented_name)),
            prettytime(trial.time, true), prettymemory(trial.memory, true))
    end

    print(divider)

    print_benchmark_details("Details", base_dict)
end

"Prints the comparison of two JSON file benchmarks to stdout."
function print_benchmark_comparison(baseline_file, alternative_file)
    base_dict = load_benchmark_from_file(baseline_file)
    base = base_dict["benchmark"][1]

    alt_dict = load_benchmark_from_file(alternative_file)
    alt = alt_dict["benchmark"][1]

    root = get_testset_hierarchy(alt)
    column_width, divider = print_table_header(root, 39, 9, 15)

    iterate_hierarchy(root) do test
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

        @printf("%s%s    %s     %s\n",
            test.indented_name,
            " "^(column_width - length(test.indented_name)),
            color(time_color(time_diff), prettytime(trial.time, true) * " (" * prettypercent(time_diff, true) * ")"),
            color(memory_color(memory_diff), prettymemory(trial.memory, true) * " (" * prettypercent(memory_diff, true) * ")"))
    end

    print(divider)

    print_benchmark_details("Baseline", base_dict)
    print_benchmark_details("Alternative", alt_dict)
end
