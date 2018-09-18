######################################################################################
# Compares two benchmark JSON files that were previously created by run_benchmark.jl
######################################################################################
push!(LOAD_PATH, @__DIR__); using TestUtils

if length(ARGS) != 2
    println("Specify two benchmark JSON files as input:\n\n$(basename(@__FILE__)) baseline_file alternative_file")
else
    compare_benchmarks(ARGS[1], ARGS[2])
end
