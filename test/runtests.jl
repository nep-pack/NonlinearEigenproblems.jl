################################################################################
# Runs all .jl files in the current directory containing a @test or @testset
# macro, except those specified to be excluded
################################################################################
using Base.Test
using TimerOutputs

# Add tests below if you wish that they are not run together with all tests
tests_not_to_run = Set{String}(map(uppercase, [
    "runtests.jl", # this file
    "Beyn_parallel.jl", # currently disabled
    "fiber.jl", # needs MATLAB
    "gun.jl", # needs MATLAB
    "matlablinsolvers.jl", # needs MATLAB
    "wep_large.jl", #  Extensive test for used during development. Needs MATLAB
    ]))

function is_test_script(file)
    if ismatch(r"(?i)\.jl$", file)
        src = open(readstring, file)
        pos = 1
        while !done(src, pos)
            expr, pos = parse(src, pos)
            if contains_test_macro(expr)
                return true
            end
        end
    end
    return false
end

function contains_test_macro(expr::Expr)
    if expr.head == :macrocall && in(expr.args[1], [Symbol("@test") Symbol("@testset")])
        return true
    end
    return any(e -> contains_test_macro(e), filter(a -> isa(a, Expr), expr.args))
end

include("load_modules_for_tests.jl")

@testset "All tests" begin
    base_path = string(@__DIR__)
    file_list = readdir(base_path)
    tests_to_run = filter(f -> is_test_script(joinpath(base_path, f)) && !in(uppercase(f), tests_not_to_run), file_list)

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
