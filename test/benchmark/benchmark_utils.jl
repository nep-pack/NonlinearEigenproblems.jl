push!(LOAD_PATH, normpath(@__DIR__, "..")); using TestUtils
using InteractiveUtils
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
    if haskey(TestUtils.benchmark_results, test)
        trial = minimum(TestUtils.benchmark_results[test])
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
function save_benchmark(test_results, file_name)
    aggregated_trials = aggregate_trials(Dict(), test_results, [])
    benchmark = BenchmarkGroup()
    foreach(t -> benchmark[t.key] = t.trial, values(aggregated_trials))

    dict = Dict()

    dict["time"] = Dates.now(Dates.UTC)
    dict["config"] = @sprintf("seconds=%.1f", TestUtils.benchmark_duration_seconds)
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
load_benchmark_file(file_name) = load_benchmark(JSON.parsefile(file_name))

"Load benchmark results from a JSON string"
load_benchmark_string(str) = load_benchmark(JSON.parse(str))

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
