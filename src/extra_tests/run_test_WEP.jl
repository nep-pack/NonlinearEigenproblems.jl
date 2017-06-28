workspace()
push!(LOAD_PATH, pwd())	# looks for modules in the current directory
using Gallery

#OBS: Only needed to run the debug
push!(LOAD_PATH, pwd()*"/gallery_extra/waveguide")	# looks for modules in the correct directory
using waveguide_debug



println("===========================")
println("||   This is WEP-test    ||")
println("===========================")

delta = 0.1

matlab_debug_eigval_comp_WEP_FD_and_SPMF(105, 35, delta)

debug_WEP_FD_preconditioner(delta)

matlab_debug_WEP_FD(119, 115, delta)

matlab_debug_full_matrix_WEP_FD_SPMF(21, 17, delta)

debug_sqrtm_schur(281)

debug_sqrt_derivative()

debug_Mlincomb_FD_WEP(29, 31, delta)

matlab_debug_Schur_WEP_FD_SPMF(49, 45, delta)

fft_debug_mateq(431, 427, delta)

debug_Sylvester_SMW_WEP(109, 105, delta, 7)
