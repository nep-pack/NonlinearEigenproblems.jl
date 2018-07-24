##################
# Runs all tests #
##################
using Base.Test
using TimerOutputs

function to_uppercase_set(array)
    return Set{String}(map(uppercase, array))
end

# Add tests below if you wish that they are not run together with all tests
tests_not_to_run = to_uppercase_set([
    "runtests.jl", # this file
    "submit_test_coverage.jl", # run by CI tool
    "Beyn_parallel.jl", # currently disabled
    "fiber.jl", # needs MATLAB
    "gun.jl", # needs MATLAB
    "matlablinsolvers.jl", # needs MATLAB
    ])

@testset "All tests" begin
    base_path = string(@__DIR__)
    file_list = readdir(base_path)
    tests_to_run = filter(f -> ismatch(r"(?i)\.jl$", f) && !in(uppercase(f), tests_not_to_run), file_list)

    to = TimerOutput()

    for i = 1:length(tests_to_run)
        file = tests_to_run[i]
        test_name = replace(file, r"(?i)\.jl$", "")
        @printf("Running test %s (%d / %d)\n", test_name, i, length(tests_to_run))
        @timeit to test_name include(base_path * "/" * file)
    end

    show(to; title = "Test Performance", compact = true)
    println()
end
