#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery

println("===========================")
println("||   This is WEP-test    ||")
println("===========================")



nep = nep_gallery("waveguide", 10, 11, "JaRLeBring", "FD", "SPMF", 0.1)

nep





