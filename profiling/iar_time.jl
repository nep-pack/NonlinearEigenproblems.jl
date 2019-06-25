#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

workspace(); push!(LOAD_PATH, string(@__DIR__, "/../src"));push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra")); push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide")); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver; using Gallery; using IterativeSolvers; using Test


#import NEPSolver.inner_solve; include("../src/inner_solver.jl");import NEPSolver.iar; include("../src/method_iar.jl");import NEPSolver.tiar; include("../src/method_tiar.jl");import NEPSolver.iar_chebyshev; include("../src/method_iar_chebyshev.jl");


nep=nep_gallery("dep0_tridiag",1000)
import NEPCore.compute_Mlincomb
function compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))
	n,k=size(V); Av=get_Av(nep)
	V=broadcast(*,V,a.');
	T=eltype(V)
	z=zeros(T,n)
	for j=1:length(nep.tauv)
		w=exp(-λ*nep.tauv[j])*(-nep.tauv[j]).^(0:k-1);
		z+=Av[j+1]*sum(broadcast(*,V,w.'),2);
	end
	if k>1 z-=view(V,:,2:2) end
	return z-λ*view(V,:,1:1)
end

iar(nep,neigs=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
@time iar(nep,neigs=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)

iar_chebyshev(nep,neigs=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
@time iar_chebyshev(nep,neigs=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)

tiar(nep,neigs=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
@time tiar(nep,neigs=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=100)
1
