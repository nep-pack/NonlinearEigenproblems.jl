#Intended to be run from nep-pack/ directory or nep-pack/profiling directory
using NonlinearEigenproblems

using LinearAlgebra
using IterativeSolvers
using Random
using SparseArrays
import Base.exp


n=100;
TT=Complex{Float16};
Random.seed!(1) # reset the random seed
K=[1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]; # sparsity pattern of tridiag matrix
A1=sparse(K, J, rand(Complex{Float32},3*n-2))
A2=sparse(K, J, rand(Complex{Float32},3*n-2))
A3=sparse(K, J, rand(Complex{Float32},3*n-2))
A4=sparse(K, J, rand(Complex{Float32},3*n-2))

AA = [A1,A2,A3,A4]
nep=PEP(AA)

n=size(nep,1);	k=5;
V=rand(TT,n,k);	λ=rand(TT)*im+rand(TT); a=rand(TT,k)

eltype_nep=eltype(nep.A[1])
predict_type=promote_type(promote_type(typeof(λ),eltype(V)),eltype_nep)


z1=compute_Mlincomb(nep,λ,V)
eltype(z1)==predict_type
