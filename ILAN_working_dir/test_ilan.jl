using NonlinearEigenproblems, Random, SparseArrays, Test, LinearAlgebra
import ..NEPSolver.ilan;
import ..NEPSolver.ilan_benchmark;
include("../src/method_ilan.jl");
include("../src/method_ilan_benchmark.jl");

n=100
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';
A4 = sparse(K, J, rand(3*n-2)); A4 = A4+A4';

f1= S -> one(S)
f2= S -> -S
f3= S -> exp(-S)
f4= S -> sqrt(one(S)-2*S)

v0=rand(n)
nep=SPMF_NEP([A1,A2,A3,A4],[f1,f2,f3,f4])

λ,v=ilan(nep,σ=0,γ=1;neigs=10,logger=1,maxit=40,tol=eps()*100,check_error_every=5,v=v0)
