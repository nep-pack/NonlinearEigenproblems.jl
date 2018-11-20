using NonlinearEigenproblems, Random, SparseArrays, Revise, Profile

import ..NEPSolver.ilan;
include("../src/method_ilan.jl");
Profile.clear()
Profile.init(n = 10^7, delay = 0.01)

n=1000
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';

nep=DEP([A1,A2,A3],[0,1,0.8])

σ=0;
γ=1;

#Dnep=DerSPMF(nep,σ,400)
ilan(nep;σ=σ,γ=γ,Neig=10,displaylevel=1,maxit=200,tol=eps()*100,check_error_every=1)
@profile ilan(nep;σ=σ,γ=γ,Neig=10,displaylevel=1,maxit=200,tol=eps()*100,check_error_every=1)
Profile.print(format=:flat,sortedby=:count)
