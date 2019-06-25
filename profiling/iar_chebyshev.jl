#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

#workspace(); push!(LOAD_PATH, string(@__DIR__, "/../src"));push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra")); push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide")); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver; using Gallery; using IterativeSolvers; using Test

import NEPSolver.inner_solve; include("../src/inner_solver.jl");import NEPSolver.iar_chebyshev; include("../src/method_iar_chebyshev.jl");

Profile.clear()
nep=nep_gallery("dep0_tridiag",5000)
Profile.init(n = 10^7, delay = 0.01)
@profile iar_chebyshev(nep,neigs=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
#Profile.print()
Profile.print(format=:flat,sortedby=:count)
