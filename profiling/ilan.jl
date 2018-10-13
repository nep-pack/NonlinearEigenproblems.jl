using NonlinearEigenproblems, Random, SparseArrays, Revise
include("../src/inner_solver.jl");
include("../src/method_ilan.jl");


n=4;
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2))
A2 = sparse(K, J, rand(3*n-2))
A3 = sparse(K, J, rand(3*n-2))
A1 = A1+A1';
A2 = A2+A2';
A3 = A3+A3';


f1= S -> one(S)
f2= S -> -S
f3= S -> exp(-S)

nep=SPMF_NEP([A1,A2,A3],[f1,f2,f3])

V=ilan(nep,σ=0,γ=1;Neig=10,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)
V
