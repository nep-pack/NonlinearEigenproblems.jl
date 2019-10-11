using NonlinearEigenproblems, Random, SparseArrays, Test, LinearAlgebra, PyPlot, Profile, BenchmarkTools
import ..NEPSolver.ilan;
import ..NEPSolver.iar;


include("../src/method_ilan.jl");


n=10000
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';
nep=DEP([A1,A2,A3],[0,1,1])
v0=rand(n)

@btime λ2,v2,err,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=1,maxit=50,tol=1e-6,check_error_every=Inf,v=v0)
1
