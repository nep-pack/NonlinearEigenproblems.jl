#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

#workspace(); push!(LOAD_PATH, string(@__DIR__, "/../src"));push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra")); push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide")); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver; using Gallery; using IterativeSolvers; using Base.Test




nep=nep_gallery("dep0_tridiag",10000)
#nep=nep_gallery("dep0")
#nep.tauv[2]=0
#nep.A[1]=zeros(nep.A[1])
#nep.A[2]=zeros(nep.A[1])



function compute_Mlincomb_DEP(nep::DEP,λ::Number,V;a=ones(size(V,2)))
	# TODO: fix this function in a way that computes the derivatives also for \lambda != 0
	# now it works for real \lambda
	#nep.Md_lin_comb=@(X,j) -X(:,1)+A1*(sum(bsxfun(@times, X(:,1:j),(-1).^(1:j)),2));
	k=size(V,2)
	Av=get_Av(nep)
	z=zeros(V[:,1])
	for j=1:length(nep.tauv)
		w=exp(-λ*nep.tauv[j])*(-nep.tauv[j]).^(0:k-1);
		z+=Av[j+1]*sum(broadcast(*,V,w.'),2);
	end
	if k>1
		z=-V[:,2]+z
	end
		z=z-λ*V[:,1]
	return z
end

#import NEPCore.compute_Mlincomb
#compute_Mlincomb(nep::DEP,λ::Number,V;a=ones(size(V,2)))=compute_Mlincomb_DEP(nep::DEP,λ::Number,V;a=ones(size(V,2)))

n=size(nep,1);	k=1;
V=rand(n,k);	λ=im;	#TODO: if λ complex doesn't work. WHY?
z1=compute_Mlincomb_DEP(nep,λ,V)
@time z1=compute_Mlincomb_DEP(nep,λ,V)


import NEPCore.compute_Mlincomb
z2=compute_Mlincomb(nep,λ,V)
@time z2=compute_Mlincomb(nep,λ,V)


println("Error=",norm(z1-z2))
