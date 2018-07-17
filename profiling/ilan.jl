#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

#workspace(); push!(LOAD_PATH, string(@__DIR__, "/../src"));push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra")); push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide")); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver; using Gallery; using IterativeSolvers; using Base.Test

import NEPSolver.inner_solve; include("../src/inner_solver.jl");import NEPSolver.ilan; include("../src/method_ilan.jl");

nep=nep_gallery("dep0",100)
ilan(nep,σ=-1,γ=2;Neig=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
