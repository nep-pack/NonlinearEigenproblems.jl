#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

workspace(); push!(LOAD_PATH, string(@__DIR__, "/../src"));push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra")); push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide")); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver; using Gallery; using IterativeSolvers; using Base.Test


#import NEPSolver.inner_solve; include("../src/inner_solver.jl");import NEPSolver.iar; include("../src/method_iar.jl");import NEPSolver.tiar; include("../src/method_tiar.jl");import NEPSolver.iar_chebyshev; include("../src/method_iar_chebyshev.jl");


#nep=nep_gallery("dep0_tridiag",10000)
nep=nep_gallery("dep0")
import NEPCore.compute_Mlincomb

function compute_Mlincomb_DEP(nep::DEP,λ::Number,V;a=ones(size(V,2)))
	return V[:,1]
end

compute_Mlincomb(nep::NEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_DEP(nep::NEP,λ::Number,V;a=ones(size(V,2)))

n=size(nep,1);	k=3;
V=rand(n,k);	λ=-1+1im;
v[:,1]=ones(size(nep,1));
compute_Mlincomb(nep,λ,v)
