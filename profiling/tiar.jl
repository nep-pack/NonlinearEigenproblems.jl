#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

#workspace(); push!(LOAD_PATH, string(@__DIR__, "/../src"));push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra")); push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide")); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver; using Gallery; using IterativeSolvers; using Test

import NEPSolver.inner_solve; include("../src/inner_solver.jl");import NEPSolver.tiar; include("../src/method_tiar.jl");

Profile.clear()
nep=nep_gallery("dep0_tridiag",1000)
import NEPCore.compute_Mlincomb
#function compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))
#    return compute_Mlincomb_from_Mder(nep,λ,V,a)
#end
Profile.init(n = 10^7, delay = 0.01)
@profile tiar(nep,σ=-1,γ=2;Neig=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
#Profile.print()
Profile.print(format=:flat,sortedby=:count)
