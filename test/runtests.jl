##################
# Runs all tests #
##################
using Base.Test

function to_uppercase_set(array)
    return Set{String}(map(uppercase, array))
end

# Add tests below if you wish that they are not run together with all tests
tests_not_to_run = to_uppercase_set([
    "runtests.jl", # this file
    "Beyn_parallel.jl", # currently disabled
    "fiber.jl", # needs MATLAB
    "gun.jl", # needs MATLAB
    "matlablinsolvers.jl", # needs MATLAB
    "iar_chebyshev.jl", # currently contains a bug
    ])

@testset "All tests" begin
    base_path = string(@__DIR__)
    file_list = readdir(base_path)
    tests_to_run = filter!(f -> ismatch(r"(?i)\.jl$", f) && !in(uppercase(f), tests_not_to_run), file_list)

    for i = 1:length(tests_to_run)
        file = tests_to_run[i]
        @printf("Running test %s (%d / %d)\n", file, i, length(tests_to_run))
        @time include(base_path *"/" * file)
    end
end
