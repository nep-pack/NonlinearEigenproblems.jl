#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

workspace(); push!(LOAD_PATH, string(@__DIR__, "/../src"));push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra")); push!(LOAD_PATH, string(@__DIR__, "/../src/gallery_extra/waveguide")); using NEPCore; using NEPTypes; using LinSolvers; using NEPSolver; using Gallery; using IterativeSolvers; using Base.Test

nep=nep_gallery("dep0_tridiag",1000)


n=size(nep,1);	k=60;
V=rand(n,k);	λ=rand()*im+rand();	#TODO: if λ complex doesn't work. WHY?
a=rand(k)

z1=compute_Mlincomb(nep,λ,copy(V),a=a)

# structure exploitation for DEP (TODO: document this)
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

@time z1=compute_Mlincomb(nep,λ,copy(V),a=a)

# old way of compute_Mlincomb used for DEP
import NEPCore.compute_Mlincomb_from_MM
z2=compute_Mlincomb_from_MM(nep,λ,V,a)
@time z2=compute_Mlincomb_from_MM(nep,λ,V,a)


println("Error=",norm(z1-z2))
