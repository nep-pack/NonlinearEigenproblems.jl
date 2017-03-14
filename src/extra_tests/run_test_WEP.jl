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


nx = 179;
nz = 131;
delta = 0.1;
nep1 = nep_gallery("waveguide", nx, nz, "TAUSCH", "fD", "SPMF", delta)



nep2 = nep_gallery("waveguide", nx, nz, "jArleBRIng", "fD", "SpmF", delta)



matlab_debug_WEP_FD(nx, nz, delta)



