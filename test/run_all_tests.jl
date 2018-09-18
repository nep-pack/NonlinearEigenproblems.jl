push!(LOAD_PATH, @__DIR__); using TestUtils
using Printf
using Test

# Add tests below if you wish that they are not run together with all tests
const tests_not_to_run = Set{String}(map(uppercase, [
    "run_all_tests.jl", # this file
    "beyn_parallel.jl", # currently disabled
    "fiber.jl", # needs MATLAB
    "gun.jl", # needs MATLAB
    "wep_large.jl", #  Extensive test for used during development. Needs MATLAB
    "nleigs_test_utils.jl", # utilities used by other tests
]))

function run_all_tests(test_name_regex = "")
    @testset "All tests" begin
        root = string(@__DIR__)
        tests_to_run = [joinpath(dir, file) for (dir, _, files) in walkdir(root) for file in files
            if occursin(Regex(test_name_regex), file) && is_test_script(joinpath(dir, file)) && !in(uppercase(file), tests_not_to_run)]

        for i = 1:length(tests_to_run)
            file = tests_to_run[i]
            test_name = replace(file, Regex("$root/?(.+).jl\$", "i") => s"\1")
            @printf("Running test %s (%d / %d)\n", test_name, i, length(tests_to_run))
            @time include(file)
        end
    end
end
