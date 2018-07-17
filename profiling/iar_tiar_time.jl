#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

workspace(); push!(LOAD_PATH, string(@__DIR__, "/../src"));push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra")); push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide")); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver; using Gallery; using IterativeSolvers; using Base.Test


#import NEPSolver.inner_solve; include("../src/inner_solver.jl");import NEPSolver.iar; include("../src/method_iar.jl");import NEPSolver.tiar; include("../src/method_tiar.jl");


nep=nep_gallery("dep0_tridiag",5000)
tiar(nep,σ=-1,γ=2;Neig=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)


@time tiar(nep,σ=-1,γ=2;Neig=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)

iar(nep,σ=-1,γ=2;Neig=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
@time iar(nep,σ=-1,γ=2;Neig=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
1
