using Printf
using Test

# Add tests below if you wish that they are not run together with all tests
const tests_not_to_run = Set{String}(map(uppercase, [
    "run_all_tests.jl", # this file
    "fiber.jl", # needs MATLAB
    "gun.jl", # needs MATLAB
    "cd_player.jl", # needs MATLAB
    "wep_large.jl", #  Extensive test for used during development. Needs MATLAB
    "dtn_dimer.jl", #  Needs additional files
    "deflation2.jl", # under development
    "runtests.jl",
    "NonlinearEigenproblemsTest.jl"
]))

function run_all_tests(test_name_regex = "")
    @testset "All tests" begin
        root = string(@__DIR__)
        tests_to_run = [joinpath(dir, file) for (dir, _, files) in walkdir(root) for file in files
            if occursin(Regex(test_name_regex), file) && is_test_script(dir, file) && !in(uppercase(file), tests_not_to_run)]


        for i = 1:length(tests_to_run)
            file = tests_to_run[i]
            test_name = replace(file, Regex("$root/?(.+).jl\$", "i") => s"\1")
            @printf("Running test %s (%d / %d)\n", test_name, i, length(tests_to_run))
            @time include(file)
        end
    end
end

function is_test_script(dir::AbstractString, file::AbstractString)
    # match files ending in ".jl" and not containing # or ~ (emacs temp files)
    if occursin(r"(?i)^[^#~]+\.jl$", file)
        src = read(joinpath(dir, file), String)

        pos = 1
        while pos <= length(src)
            expr, pos = Meta.parse(src, pos)
            expr === nothing && continue
            if expr.head == :incomplete
                msg = join(map(string, expr.args), "; ")
                throw(Meta.ParseError("While parsing file $file got invalid expression: $msg"))
            end
            contains_test_macro(expr) && return true
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
