using NonlinearEigenproblems, Random, SparseArrays, Revise, PyPlot
import ..NEPSolver.ilan;
include("../src/method_ilan.jl");

n=1000
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A0 = sparse(K, J, rand(3*n-2));
A1 = sparse(K, J, rand(3*n-2));

nep=DEP([A0,A1],[0,1])
λ,v=iar(nep;tol=1e-12,neigs=10);
