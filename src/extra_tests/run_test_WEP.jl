#  This is the first code in NEP-pack
workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using NEPSolver
using NEPCore
using NEPTypes
using Gallery
using LinSolvers

#OBS: Only needed to run the debug
push!(LOAD_PATH, pwd()*"/../gallery_extra/waveguide")	# looks for modules in the correct directory
using Waveguide



println("===========================")
println("||   This is WEP-test    ||")
println("===========================")



##################################################################################################################
nz = 105
nx = nz + 4
delta = 0.1
N = 35

println("\n     Test Tausch FD WEP-format\n")
nep_jar_wep = nep_gallery("waveguide", nx, nz, "Jarlebring", "fD", "weP", delta)


σ = -0.5-0.4im
println("Generating preconditioner")
precond = @time generate_preconditioner(nep_jar_wep, N, σ)

#gmres_kwargs = ((:maxiter,100), (:restart,100), (:log,false), (:Pl,precond), (:verbose,true), (:tol, 1e-13))
gmres_kwargs = ((:maxiter,100), (:restart,100), (:log,false), (:Pl,precond), (:tol, 1e-13))
function wep_gmres_linsolvercreator(nep::NEP, λ)
    return wep_linsolvercreator(nep, λ, gmres_kwargs)
end


λ_jar_wep=NaN
x_jar_wep=NaN
try
    λ_jar_wep,x_jar_wep =resinv(nep_jar_wep, displaylevel=1, λ=σ, maxit = 25, tolerance = 1e-10, v=ones(Complex128,nx*nz+2*nz), c=zeros(Complex128,nx*nz+2*nz), linsolvercreator=wep_gmres_linsolvercreator)
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ_jar_wep = e.λ
    x_jar_wep = e.v
end
println("Resnorm: ", compute_resnorm(nep_jar_wep, λ_jar_wep, x_jar_wep))
println("Eigenvalue: ", λ_jar_wep)
println("Eigenvector norm: ", norm(x_jar_wep), "\n")



##################################################################################################################
nz = 105
nx = nz + 4
delta = 0.1

println("\n     Test Tausch FD SPFM\n")
nep_tausch = nep_gallery("waveguide", nx, nz, "TAUSCH", "fD", "SPMF", delta)

λ_tausch=NaN
x_tausch=NaN
try
    λ_tausch,x_tausch =resinv(nep_tausch;displaylevel=1, λ=-0.015-5.1im, maxit = 25, tolerance = 1e-10, v=ones(Complex128,nx*nz+2*nz), c=zeros(Complex128,nx*nz+2*nz))
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ_tausch=e.λ
    x_tausch=e.v
end
println("Resnorm: ", compute_resnorm(nep_tausch,λ_tausch,x_tausch))
println("Eigenvalue: ", λ_tausch)
println("Eigenvector norm: ", norm(x_tausch), "\n")


##################################################################################################################
nz = 105
nx = nz + 4
delta = 0.1

println("\n     Test Jarlebring FD SPFM\n")
nep_jar = nep_gallery("waveguide", nx, nz, "jArleBRIng", "fD", "SpmF", delta)

λ_jar=NaN
x_jar=NaN
try
    λ_jar,x_jar =resinv(nep_jar;displaylevel=1, λ=-0.5-0.4im, maxit = 25, tolerance = 1e-10, v=ones(Complex128,nx*nz+2*nz), c=zeros(Complex128,nx*nz+2*nz))
catch e
    # Only catch NoConvergence
    isa(e, NoConvergenceException) || rethrow(e)
    println("No convergence because:"*e.msg)
    # still access the approximations
    λ_jar=e.λ
    x_jar=e.v
end
println("Resnorm: ", compute_resnorm(nep_jar,λ_jar,x_jar))
println("Eigenvalue: ", λ_jar)
println("Eigenvector norm: ", norm(x_jar), "\n")


##################################################################################################################
delta = 0.1

println("\n     Different DEBUG-tests\n")

matlab_debug_WEP_FD(119, 115, delta)

matlab_debug_full_matrix_WEP_FD_SPMF(21, 17, delta)

debug_sqrtm_schur(281)

debug_sqrt_derivative()

debug_Mlincomb_FD_WEP(29, 31, delta)

matlab_debug_Schur_WEP_FD_SPMF(49, 45, delta)

fft_debug_mateq(431, 427, delta)

debug_Sylvester_SMW_WEP(109, 105, delta, 7)

