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
f4 = S -> sqrt(10*S.+100*I);

n=100000;
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
V=rand(n,k);	λ=rand()*im+rand();
a=rand(k)

z1=compute_Mlincomb(nep,λ,copy(V),a)
@time z1=compute_Mlincomb(nep,λ,V,a)
# old way of compute_Mlincomb used for DEP
import NonlinearEigenproblems.NEPCore.compute_Mlincomb_from_MM
z2=compute_Mlincomb_from_MM(nep,λ,V,a)
@time z2=compute_Mlincomb_from_MM(nep,λ,V,a)
println("Error=",norm(z1-z2))
