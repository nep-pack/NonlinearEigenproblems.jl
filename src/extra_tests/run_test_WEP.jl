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


nx = 137;
nz = 11;
delta = 0.1;
nep_tausch = nep_gallery("waveguide", nx, nz, "TAUSCH", "fD", "SPMF", delta)

#λ_tausch=NaN;
#x_tausch=NaN
#try
#    λ_tausch,x_tausch =newton(nep_tausch,displaylevel=1, λ=-10-10im, maxit = 50, tolerance = 1e-10);
#catch e
#    # Only catch NoConvergence
#    isa(e, NoConvergenceException) || rethrow(e)
#    println("No convergence because:"*e.msg)
#    # still access the approximations
#    λ_tausch=e.λ
#    x_tausch=e.v
#end
#println("Resnorm: ",compute_resnorm(nep_tausch,λ_tausch,x_tausch))
#println("Eigenvalue: ", λ_tausch)
#println("Eigenvector norm: ", norm(x_tausch))


nep_jar = nep_gallery("waveguide", nx, nz, "jArleBRIng", "fD", "SpmF", delta)

#λ_jar=NaN;
#x_jar=NaN
#try
#    λ_jar,x_jar =newton(nep_jar,displaylevel=1, λ=-10-10im, maxit = 50, tolerance = 1e-10);
#catch e
#    # Only catch NoConvergence
#    isa(e, NoConvergenceException) || rethrow(e)
#    println("No convergence because:"*e.msg)
#    # still access the approximations
#    λ_jar=e.λ
#    x_jar=e.v
#end
#println("Resnorm: ",compute_resnorm(nep_jar,λ_jar,x_jar))
#println("Eigenvalue: ", λ_jar)
#println("Eigenvector norm: ", norm(x_jar))

debug_sqrtm_schur(281)

matlab_debug_WEP_FD(nx, nz, delta)



