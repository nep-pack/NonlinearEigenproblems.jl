using BenchmarkTools
using TimeZones
using Dates
using HTTP
using JSON

include("benchmark_utils.jl")

const CI_GITHUB_REPO = "maxbennedich/julia-ci"
const SOURCE_GITHUB_REPO = "nep-pack/NonlinearEigenproblems.jl"
const NR_COMMITS_TO_SHOW = 16

abstract type BenchmarkFile end

struct GitHubBenchmarkFile <: BenchmarkFile
    build_nr
    url
end

struct LocalBenchmarkFile <: BenchmarkFile
    build_nr
    path
end

build_nr(file_name) = parse(Int, split(file_name, "-")[end-1])

GitHubBenchmarkFile(entry) =
    GitHubBenchmarkFile(build_nr(entry["name"]), entry["download_url"])

LocalBenchmarkFile(file_name) =
    LocalBenchmarkFile(build_nr(file_name), file_name)

get_benchmark_json(file::GitHubBenchmarkFile) = (println("Loading $file"); String(HTTP.request("GET", file.url).body))
get_benchmark_json(file::LocalBenchmarkFile) = read(file.path, String)

import Base.show

show(io::IO, f::GitHubBenchmarkFile) = print(io, f.url)
show(io::IO, f::LocalBenchmarkFile) = print(io, f.path)

function get_github_files(branch)
    req = HTTP.request("GET", "https://api.github.com/repos/$CI_GITHUB_REPO/contents", [("User-Agent", "julia-ci")])
    json = JSON.parse(String(req.body))
    return [GitHubBenchmarkFile(f) for f in json if f["type"] == "file" && startswith(f["name"], "benchmark-$branch-")]
end

function get_local_files(branch)
    return [LocalBenchmarkFile(f) for f in readdir() if !isdir(f) && startswith(f, "benchmark-$branch-")]
end

function print_report(io, branch)
    print_report_header(io)

#    files = get_github_files(branch)
    files = get_local_files(branch)
    sort!(files, by = x -> x.build_nr)
    files = files[max(1, end-NR_COMMITS_TO_SHOW+1):end]
    println("Loaded $(length(files)) files")

    # table header
    println(io, "<table>")
    println(io, "<tr>")
    println(io, "<th class=\"test-name\">Test</th>")

    println(io, "<th></th>")
    println(io, "<th colspan=\"2\">Reference</th>")

    runs = map(f -> load_benchmark_string(get_benchmark_json(f)), files)
    bruns = map(r -> r["benchmark"][1], runs)
    println("Loaded $(length(runs)) runs")

    for run in runs
        println(io, "<th></th>")
        commit_id = split(run["git"], " ")[1]
        gitstr = run["git"]
        length(gitstr) > 160 && (gitstr = gitstr[1:157] * "...")
        tooltip = "<b>Run:</b> $(run["time"]); $(run["config"])<br/>" *
            "<b>CPU:</b> $(run["cpu"])<br/>" *
            "<b>Git:</b> $gitstr"
        id = "<a href=\"https://github.com/$SOURCE_GITHUB_REPO/commit/$commit_id\">$commit_id</a>";
        println(io, "<th colspan=\"2\"><div class=\"tooltip\"><span class=\"tooltiptext\">$tooltip</span>$id</div></th>")
    end

    println(io, "</tr>")

    root = get_testset_hierarchy(bruns[end])

    iterate_hierarchy(root) do test
        println(io, "<tr>")
        spaces = findfirst(!isequal(' '), test.indented_name)
        indented_name = "&nbsp;"^((spaces-1)*2) * test.indented_name[spaces:end]
        println(io, "<td class=\"test-name\">$indented_name</td>")

        key = test.trial[1]

        min_time = mapreduce(run -> haskey(run, key) ? run[key].time : Inf, min, bruns)
        min_memory = mapreduce(run -> haskey(run, key) ? run[key].memory : Inf, min, bruns)

        println(io, "<td></td>")
        println(io, "<td>$(prettytime(min_time))</td>")
        println(io, "<td>$(prettymemory(min_memory))</td>")

        for run in bruns
            println(io, "<td></td>")

            if haskey(run, key)
                b = run[key]

                time_increase_pct = 100 * (b.time / min_time - 1)
                time_str = prettypercent(time_increase_pct)
                time_class = time_increase_pct >= 10 ? "bad" : time_increase_pct >= 5 ? "warn" : "ok"

                memory_increase_pct = 100 * (b.memory / min_memory - 1)
                memory_str = prettypercent(memory_increase_pct)
                memory_class = memory_increase_pct >= 5 ? "bad" : memory_increase_pct >= 1 ? "warn" : "ok"
            else
                time_str = prettypercent(NaN)
                memory_str = prettypercent(NaN)
                time_class = "missing"
                memory_class = "missing"
            end

            println(io, "<td class=\"$time_class\">$time_str</td>")
            println(io, "<td class=\"$memory_class\">$memory_str</td>")
        end

        println(io, "</tr>")
    end

    println(io, "</table>")

    print_report_footer(io)
end

function print_report_header(io)
    println(io, """
<!DOCTYPE html>
<html>
<head>
<meta charset='utf-8'>
<title>Regression Analysis</title>
<style>
@import url(http://fonts.googleapis.com/css?family=Roboto+Condensed:400,700);
html * {
  font-family: 'Roboto Condensed', sans-serif;
  margin:0;
  padding:0;
  font-size:15px;
}
a {
  color: #428bca;
  text-decoration: none;
}
a:hover, a:focus {
  color: #2a6496;
  text-decoration: underline;
}
#header {
  background-color: #2f2f2f;
  width: 100%;
  padding-bottom: 20px;
  padding-top: 20px;
  color:#fff;
  text-align:center;
  font-weight: 700;
  font-size: 48px;
}
table {
  display: inline-block;
  overflow-x: auto;
  white-space: nowrap;
}
th, td {
  padding-left: 5px;
  padding-right: 5px;
  padding-top: 0px;
  padding-bottom: 0px;
  text-align: center;
}
.bad {
  background-color: #f99;
}
.warn {
  background-color: #fd7;
}
.ok {
  background-color: #9d9;
}
.missing {
  background-color: #fff;
}
.test-name {
  text-align: left;
}
.tooltip {
  position: relative;
  display: inline-block;
}
.tooltip .tooltiptext {
  font-weight: 400;
  text-align: left;

  visibility: hidden;
  top: 100%;
  right: 0%;
  background-color: rgba(0,0,0,0.75);
  color: #fff;
  border-radius: 6px;
  padding: 10px;

  position: absolute;
  z-index: 1;
}
.tooltip:hover .tooltiptext {
  visibility: visible;
}
</style>
</head>
<body>
<div id="header">Regression Analysis</div><br/>""")
end

function print_report_footer(io)
    println(io, """
</body>
</html>""")
end

open("test/benchmark_report.html", "w") do io
    print_report(io, "benchmarks")
end
