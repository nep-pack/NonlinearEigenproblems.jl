################################################################################
# Runs the test suite in benchmark mode and saves the results to a JSON file
################################################################################
using NonlinearEigenproblemsTest
using Logging

# turn off logging of @info statements and test output
Base.CoreLogging.disable_logging(Logging.Info)
set_displaylevel(0)

include("benchmark_utils.jl")
include(normpath("..", "run_all_tests.jl"))

if length(ARGS) < 1
    println("Benchmark JSON file name required:\n\n$(basename(@__FILE__)) benchmark_file [duration_seconds] [test_name_regex]")
else
    enable_benchmark()
    length(ARGS) >= 2 && set_benchmark_duration_seconds(parse(Float64, ARGS[2]))
    test_name_regex = length(ARGS) < 3 ? "" : ARGS[3]
    test_results = @timed run_all_tests(test_name_regex)
    save_benchmark(test_results[1:2]..., ARGS[1])
    print_benchmark(ARGS[1])
end
