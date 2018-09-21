################################################################################
# Runs the test suite in benchmark mode and saves the results to a JSON file
################################################################################
push!(LOAD_PATH, @__DIR__); using TestUtils
using Logging

# turn off logging of @info statements
Base.CoreLogging.disable_logging(Logging.Info)

include("run_all_tests.jl")

if length(ARGS) < 1
    println("Benchmark JSON file name required:\n\n$(basename(@__FILE__)) benchmark_file [duration_seconds] [test_name_regex]")
else
    enable_benchmark()
    length(ARGS) >= 2 && set_benchmark_duration_seconds(parse(Float64, ARGS[2]))
    test_name_regex = length(ARGS) < 3 ? "" : ARGS[3]
    save_benchmark(run_all_tests(test_name_regex), ARGS[1])
end
