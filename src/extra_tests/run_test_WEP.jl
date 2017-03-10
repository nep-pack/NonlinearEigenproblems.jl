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


nx = 10;
nz = 11;
delta = 0.1;
nep1 = nep_gallery("waveguide", nx, nz, "TAUSCH", "fD", "SPMF", delta)

using MATLAB
WEP_path = "../matlab/WEP"
@mput nx nz delta WEP_path
@matlab begin
    addpath(WEP_path)
    nxx = double(nx)
    nzz = double(nz)
    options = struct
    options.delta = delta
    options.wg = "TAUSCH"
    matlab_nep = nep_wg_generator(nxx, nzz, options)
    
    C1_t = matlab_nep.C1
    C2T_t = matlab_nep.C2T
    K_t = matlab_nep.K


    options.wg = "CHALLENGE"
    matlab_nep = nep_wg_generator(nxx, nzz, options)

    C1_j = matlab_nep.C1
    C2T_j = matlab_nep.C2T
    K_j = matlab_nep.K

    @matlab end
@mget K_t C2T_t C1_t K_j C2T_j C1_j

println(norm(K_t-nep1))


nep2 = nep_gallery("waveguide", nx, nz, "jArleBRIng", "fD", "SpmF", delta)


println(norm(K_j-nep2))




