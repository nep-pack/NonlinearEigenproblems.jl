using NonlinearEigenproblems, Random, SparseArrays, Revise
include("../src/inner_solver.jl");
include("../src/method_ilan.jl");
#include("../src/method_iar.jl");


n=200
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';
A4 = sparse(K, J, rand(3*n-2)); A4 = A4+A4';

f1= S -> one(S)
f2= S -> -S
f3= S -> exp(-S)
#f4= S -> sqrt(one(S)-2*S)
f4= S -> exp(-S)


nep=SPMF_NEP([A1,A2,A3,A4],[f1,f2,f3,f4])
V,H,ω=ilan(nep,σ=0,γ=1;Neig=10,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)
Q,_,_=svd(V)

# project (hardcoded for now)
AA1=Q'*A1*Q;
AA2=Q'*A2*Q;
AA3=Q'*A3*Q;
AA4=Q'*A4*Q;
err_lifted=(λ,z)->compute_resnorm(nep,λ,Q*z)/n;
pnep=SPMF_NEP([AA1,AA2,AA3,AA4],[f1,f2,f3,f4])
λ,_,err=iar(pnep,σ=0,γ=1;Neig=10,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1,errmeasure=err_lifted)
