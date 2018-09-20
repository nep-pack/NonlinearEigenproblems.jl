#Intended to be run from nep-pack/ directory or nep-pack/profiling directory
using NonlinearEigenproblems

using LinearAlgebra
using IterativeSolvers
using Random
using SparseArrays
import Base.exp


n=100;
Random.seed!(1) # reset the random seed
K=[1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]; # sparsity pattern of tridiag matrix
A1=sparse(K, J, rand(Float32,3*n-2))
A2=sparse(K, J, rand(Float32,3*n-2))
A3=sparse(K, J, rand(Float32,3*n-2))
A4=sparse(K, J, rand(Float32,3*n-2))

AA = [A1,A2,A3,A4]
nep=PEP(AA)

n=size(nep,1);	k=5;
V=rand(n,k);	λ=rand()*im+rand(); a=rand(k)
λ=0;



z1=compute_Mlincomb(nep,λ,copy(V),a)
@time z1=compute_Mlincomb(nep,λ,V,a)
# old way of compute_Mlincomb used for SPMF
import NonlinearEigenproblems.NEPCore.compute_Mlincomb_from_MM
z2=compute_Mlincomb_from_MM(nep,λ,V,a)
@time z2=compute_Mlincomb_from_MM(nep,λ,V,a)
println("Error=",norm(z1-z2))
