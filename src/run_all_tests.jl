#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery

println("")
println("")
println("=================================================")
println("=================================================")
println("=================================================")
println("||           THIS RUNS ALL TESTS!              ||")
println("||      Files on format 'run_test_*.jl'        ||")
println("=================================================")
println("=================================================")
println("=================================================")
println("")
println("")

sleep(2.5)



for k = 1:2
    if( k == 1)
        path = pwd()
    elseif( k == 2)
        path = path * "extra_tests/"
    end

    file_list = readdir(path)
    for i = 1:length(file_list)
        file = file_list[i]
        upp_file = uppercase(file)
        if startswith(upp_file, "RUN_TEST_") && endswith(upp_file, ".JL")
            println("=================================================")
            println("=================================================")
            println("=================================================")
            println("    Running: '" , file ,  "'")
            println("=================================================")
            println("=================================================")
            println("=================================================")
            include(path *"/" * file)
            println("")
            println("")
        end
    end
end


