#Intended to be run from nep-pack/ directory or nep-pack/profiling directory

using NonlinearEigenproblems.NEPSolver
using NonlinearEigenproblems.Gallery
using NonlinearEigenproblems.NEPCore
using Test
using LinearAlgebra
using IterativeSolvers

nep=nep_gallery("dep0_tridiag",5000000)


n=size(nep,1);	k=1;
V=rand(n,k);	λ=rand()*im+rand();
a=rand(k)

z1=compute_Mlincomb!(nep,λ,copy(V),a)

compute_Mlincomb(nep,λ,V,a)
@time z1=compute_Mlincomb!(nep,λ,V,a)
1
