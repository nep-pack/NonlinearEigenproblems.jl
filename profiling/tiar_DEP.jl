#Intended to be run from nep-pack/ directory or nep-pack/profiling directory
#workspace(); push!(LOAD_PATH, string(@__DIR__, "/../src"));push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra")); push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide")); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver; using Gallery; using IterativeSolvers; using Test


import NEPSolver.inner_solve; include("../src/inner_solver.jl");import NEPSolver.tiar; include("../src/method_tiar.jl"); import NEPSolver.iar; include("../src/method_iar.jl");

nep=nep_gallery("dep0_tridiag",1000)
#nep=nep_gallery("dep0",500)

v0=ones(size(nep,1))

tiar(nep,neigs=10,v=v0,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
@time tiar(nep,neigs=10,v=v0,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)


import NEPCore.compute_Mlincomb
compute_Mlincomb_from_MMcompute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=
compute_Mlincomb_from_MM(nep,λ,V,a)

@profile tiar(nep,neigs=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
Profile.print(format=:flat,sortedby=:count)
1
