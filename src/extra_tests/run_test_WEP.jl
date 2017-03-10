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



nep1 = nep_gallery("waveguide", 10, 11, "TAUSCH", "fD", "SPMF", 0.1)


nep2 = nep_gallery("waveguide", 10, 11, "jArleBRIng", "fD", "SpmF", 0.1)






