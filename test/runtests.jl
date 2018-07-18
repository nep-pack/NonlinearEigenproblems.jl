##################
# Run's all test #
##################
using Base.Test

# Add tests below if you wish that they are not run together with all tests
tests_not_to_run = [
    "runtests.jl", # this file
    "Beyn_parallel.jl", # currently disabled
    "fiber.jl", # needs MATLAB
    "gun.jl", # needs MATLAB
    "matlablinsolvers.jl", # needs MATLAB
    "iar_chebyshev.jl", # currently contains a bug 
    ]::Array{String,1}


base_path = string(@__DIR__)
file_list = readdir(base_path)

@testset "All tests" begin
    for i = 1:length(file_list)
        file = file_list[i]
        upp_file = uppercase(file) #OBS: Make no case difference
        if endswith(upp_file, ".JL")
            is_in_norun_list = false
            for k = 1:length(tests_not_to_run)
                is_in_norun_list = is_in_norun_list || (upp_file == uppercase(tests_not_to_run[k])) #OBS: Make no case difference
            end
            if(!is_in_norun_list)
              include(base_path *"/" * file)
            end
        end
    end
end
