#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

#workspace(); push!(LOAD_PATH, string(@__DIR__, "/../src"));push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra")); push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide")); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver; using Gallery; using IterativeSolvers; using Base.Test


import NEPSolver.inner_solve; include("../src/inner_solver.jl");import NEPSolver.tiar; include("../src/method_tiar.jl"); import NEPSolver.iar; include("../src/method_iar.jl");


nep=nep_gallery("dep0_tridiag",10000)


tiar(nep,Neig=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
@time tiar(nep,Neig=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)

import NEPCore.compute_Mlincomb
function compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))
	k=size(V,2); Av=get_Av(nep)
	V=broadcast(*,V,a.');
	z=zeros(V[:,1])
	for j=1:length(nep.tauv)
		w=exp(-λ*nep.tauv[j])*(-nep.tauv[j]).^(0:k-1);
		z+=Av[j+1]*sum(broadcast(*,V,w.'),2);
	end
	if k>1 z=-V[:,2]+z end
	return z-λ*V[:,1]
end


v0=ones(size(nep,1))
tiar(nep,Neig=10,v=v0,displaylevel=0,maxit=4,tol=eps()*100,check_error_every=100)
@time tiar(nep,Neig=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
1
