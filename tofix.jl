using NonlinearEigenproblems, Random, SparseArrays, LinearAlgebra, Revise
import ..NEPSolver.ilan;
include("src/method_ilan.jl");


n=100
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';
nep=DEP([A1,A2,A3],[0,1,0.8])
v0=rand(n)
λ,W=ilan(nep,σ=0,γ=1;neigs=Inf,logger=1,maxit=50,tol=eps()*100,check_error_every=1,v=v0,errmeasure=ResidualErrmeasure(nep),inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=50))
