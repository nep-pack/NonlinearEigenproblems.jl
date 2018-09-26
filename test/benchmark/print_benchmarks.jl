######################################################################################
# Prints a single JSON file benchmark, or the comparison of two benchmarks, that were
# previously created by run_benchmark.jl
######################################################################################
include("benchmark_utils.jl")

if !(1 <= length(ARGS) <= 2)
    println("Specify one or two benchmark JSON files as input:\n\n$(basename(@__FILE__)) baseline_file [alternative_file]")
else
    length(ARGS) > 1 ? print_benchmark_comparison(ARGS[1], ARGS[2]) : print_benchmark(ARGS[1])
end
