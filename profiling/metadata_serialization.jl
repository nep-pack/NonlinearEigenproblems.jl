using Dates
using JSON

# The below is largely copied from test/benchmark/benchmark_utils.jl

"Save metadata to a JSON file on disk"
function save_metadata(file_name, total_seconds, profile_seconds)
    dict = Dict()

    dict["time"] = Dates.now(Dates.UTC)

    dict["total-time"] = total_seconds
    dict["profile-time"] = profile_seconds

    dict["cpu-model"] = Sys.cpu_info()[1].model
    dict["cpu-name"] = Sys.CPU_NAME
    dict["cpu-threads"] = Sys.CPU_THREADS
    dict["total-memory"] = Sys.total_memory()

    dict["julia-threads"] = Threads.nthreads()

    dict["julia-version"] = string(VERSION)
    dict["julia-commit"] = Base.GIT_VERSION_INFO.commit_short

    git_commit = safe_read_command(`git log -1 --format="%h"`)
    git_branch = safe_read_command(`git rev-parse --abbrev-ref HEAD`)
    git_branch == "HEAD" && (git_branch = safe_read_command(pipeline(pipeline(`git describe --all --exact-match`, stderr=devnull), `sed 's=.*/=='`)))
    isempty(git_branch) && (git_branch = "N/A")
    git_details = safe_read_command(`git log -1 --format="%ai, %s"`)
    dict["git"] = @sprintf("%s (%s), %s", git_commit, git_branch, git_details)

    open(file_name, "w") do f
        JSON.print(f, dict)
    end

#    println("Metadata saved to file $(abspath(file_name))")
end

"Load metadata from a JSON file on disk"
function load_metadata(file_name)
    dict = JSON.parsefile(file_name)

    # format time
    dt = DateTime(dict["time"])
    dict["time"] = Dates.format(dt, "yyyy-mm-dd HH:MM:SS") * " (UTC)"

    return dict
end

function safe_read_command(cmd)
    try
        chomp(read(cmd, String))
    catch
        ""
    end
end
