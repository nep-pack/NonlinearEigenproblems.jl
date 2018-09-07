#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

using NonlinearEigenproblems.NEPSolver
using NonlinearEigenproblems.Gallery
using NonlinearEigenproblems.NEPCore
using NonlinearEigenproblems.NEPTypes
using Test
using LinearAlgebra
using IterativeSolvers

f1 = S -> -S;
f2 = S -> Matrix{eltype(S)}(I, size(S));
f3 = S -> exp(-Matrix(S));
f4 = S -> sqrt(10*S+100);

n=100;
Random.seed!(1) # reset the random seed
K=[1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]; # sparsity pattern of tridiag matrix
A1=sparse(K, J, rand(3*n-2))
A2=sparse(K, J, rand(3*n-2))
A3=sparse(K, J, rand(3*n-2))
A4=sparse(K, J, rand(3*n-2))

AA = [A1,A2,A3,A4]
fi = [f1,f2,f3,f4]
nep=SPMF_NEP(AA, fi)

n=size(nep,1);	k=50;
V=rand(n,k);	λ=rand()*im+rand();	#TODO: if λ complex doesn't work. WHY?
a=rand(k)


function compute_Mlincomb2(nep::SPMF_NEP, λ, V, a)

	n,k=size(V);
	# we need to assume that the elements of a are different than zero.
	V[:,findall(x->x==0,a)] .= 0
	a[findall(x->x==0,a)] .= 1
	S=diagm(0 => λ*ones(eltype(V),k)) + diagm(1 => (a[2:k]./a[1:k-1]).*(1:k-1))
	S=copy(transpose(S))

	Z=zeros(eltype(V),n)
	for i=1:size(nep.A,1)
		Fi=nep.fi[i](S)[:,2]
		Z=Z .+ nep.A[i]*(V*Fi);
	end
	return Z
end


z1=compute_Mlincomb2(nep,λ,copy(V),a)
compute_Mlincomb2(nep,λ,V,a)
@time z1=compute_Mlincomb2(nep,λ,V,a)

# old way of compute_Mlincomb used for DEP
import NonlinearEigenproblems.NEPCore.compute_Mlincomb_from_MM
z2=compute_Mlincomb_from_MM(nep,λ,V,a)
@time z2=compute_Mlincomb_from_MM(nep,λ,V,a)

println("Error=",norm(z1-z2))
