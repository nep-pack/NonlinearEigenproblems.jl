################################################################################
# Runs all .jl files in the directory of this script containing a @test or
# @testset macro, except those specified to be excluded in run_all_tests.jl
################################################################################
include("run_all_tests.jl")

test_name_regex = length(ARGS) < 1 ? "" : ARGS[1]
run_all_tests(test_name_regex)
