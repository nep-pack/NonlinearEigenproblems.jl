#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

using NonlinearEigenproblems.NEPSolver
using NonlinearEigenproblems.Gallery
using NonlinearEigenproblems.NEPCore
using Test
using LinearAlgebra
using IterativeSolvers

nep=nep_gallery("dep0_tridiag",10000)


n=size(nep,1);	k=1;
V=rand(n,k);	λ=rand()*im+rand();	#TODO: if λ complex doesn't work. WHY?
a=rand(k)

z1=compute_Mlincomb(nep,λ,copy(V),a)

compute_Mlincomb(nep,λ,V,a)
@time z1=compute_Mlincomb(nep,λ,V,a)

# old way of compute_Mlincomb used for DEP
import NonlinearEigenproblems.NEPCore.compute_Mlincomb_from_MM
z2=compute_Mlincomb_from_MM(nep,λ,V,a)
@time z2=compute_Mlincomb_from_MM(nep,λ,V,a)

println("Error=",norm(z1-z2))
