using NonlinearEigenproblems, Random, SparseArrays, Revise, LinearAlgebra, BenchmarkTools

n=1000
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';
A4 = sparse(K, J, rand(3*n-2)); A4 = A4+A4';

f1= S -> one(S)
f2= S -> -S
f3= S -> exp(-S)
f4= S -> exp(-S)

nep=SPMF_NEP([A1,A2,A3,A4],[f1,f2,f3,f4])
σ=0
Dnep=DerSPMF(nep,σ,200)
function time_ilan()
@btime begin ilan(Dnep;neigs=10,displaylevel=0,maxit=200,tol=eps()*100,check_error_every=1) end
@btime begin ilan(nep;neigs=10,displaylevel=0,maxit=200,tol=eps()*100,check_error_every=1) end
@btime begin tiar(nep;neigs=10,displaylevel=0,maxit=200,tol=eps()*100,check_error_every=200) end
end

time_ilan()
1
