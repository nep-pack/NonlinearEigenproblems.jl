#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

#workspace(); push!(LOAD_PATH, string(@__DIR__, "/../src"));push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra")); push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide")); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver; using Gallery; using IterativeSolvers; using Base.Test


import NEPSolver.inner_solve; include("../src/inner_solver.jl");import NEPSolver.iar; include("../src/method_iar.jl");import NEPSolver.tiar; include("../src/method_tiar.jl");import NEPSolver.iar_chebyshev; include("../src/method_iar_chebyshev.jl");


nep=nep_gallery("dep0_tridiag",100000)
#nep=nep_gallery("dep0")

function compute_Mlincomb_DEP(nep::DEP,λ::Number,V;a=ones(size(V,2)))
	#nep.Md_lin_comb=@(X,j) -X(:,1)+A1*(sum(bsxfun(@times, X(:,1:j),(-1).^(1:j)),2));
	k=size(V,2)
	Av=get_Av(nep)
	z=zeros(V[:,1])
	for j=1:length(nep.tauv)
		z+=Av[j+1]*sum(broadcast(*,V,(-nep.tauv[j]).^(1:k)'),2);
	end
	return -V[:,1]+z
end

#import NEPCore.compute_Mlincomb
#compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_DEP(nep::DEP,λ::Number,V;a=ones(size(V,2)))

n=size(nep,1);	k=100;
V=rand(n,k);	λ=-1+1im;
V[:,1]=ones(size(nep,1));
z1=compute_Mlincomb_DEP(nep,λ,V)
@time z1=compute_Mlincomb_DEP(nep,λ,V)


import NEPCore.compute_Mlincomb
z2=compute_Mlincomb_DEP(nep,λ,V)
@time z2=compute_Mlincomb_DEP(nep,λ,V)


norm(z1-z2)
