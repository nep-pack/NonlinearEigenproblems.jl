################################################################################
# Runs all .jl files in the directory of this script containing a @test or
# @testset macro, except those specified to be excluded
################################################################################
using Test
using TimerOutputs
using Printf

# Add tests below if you wish that they are not run together with all tests
tests_not_to_run = Set{String}(map(uppercase, [
    "runtests.jl", # this file
    "beyn_parallel.jl", # currently disabled
    "fiber.jl", # needs MATLAB
    "gun.jl", # needs MATLAB
    "matlablinsolvers.jl", # needs MATLAB
    "wep_large.jl", #  Extensive test for used during development. Needs MATLAB
    "nleigs_test_utils.jl", # utilities used by other tests
    ]))

function is_test_script(file::AbstractString)
    if occursin(r"(?i)\.jl$", file)
        src = read(file, String)

        pos = 1
        while pos <= length(src)
            expr, pos = Meta.parse(src, pos)
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

@testset "All tests" begin
    root = string(@__DIR__)
    tests_to_run = [joinpath(dir, file) for (dir, _, files) in walkdir(root) for file in files
        if is_test_script(joinpath(dir, file)) && !in(uppercase(file), tests_not_to_run)]

    to = TimerOutput()

    for i = 1:length(tests_to_run)
        file = tests_to_run[i]
        test_name = replace(file, Regex("$root/?(.+).jl\$", "i") => s"\1")
        @printf("Running test %s (%d / %d)\n", test_name, i, length(tests_to_run))
        @timeit to test_name include(file)
    end

    show(to; title = "Test Performance", compact = true)
    println()
end
