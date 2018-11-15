using NonlinearEigenproblems, Random, SparseArrays, Revise, LinearAlgebra, BenchmarkTools

A0=[1 3; 4 5]; A1=[3 4; 5 6];
id_op=S -> one(S) # Note: We use one(S) to be valid both for matrices and scalars
exp_op=S -> exp(S)
nep=SPMF_NEP([A0,A1],[id_op,exp_op]);
m=5
σ1=3.5;
Dnep=DerSPMF(nep,σ1,m)
σ2=9;
DDnep=DerSPMF(Dnep,σ2,m)
V=randn(2,2);
n1=compute_Mlincomb(nep,0,V)-compute_Mlincomb(Dnep,0,V)
n2=compute_Mlincomb(nep,σ1,V)-compute_Mlincomb(Dnep,σ1,V)
n3=compute_Mlincomb(nep,σ2,V)-compute_Mlincomb(Dnep,σ2,V)
n4=compute_Mlincomb(nep,0,V)-compute_Mlincomb(DDnep,0,V)
n5=compute_Mlincomb(nep,σ1,V)-compute_Mlincomb(DDnep,σ1,V)
n6=compute_Mlincomb(nep,σ2,V)-compute_Mlincomb(DDnep,σ2,V)

println("err=",norm(n1))
println("err=",norm(n2))
println("err=",norm(n3))
println("err=",norm(n4))
println("err=",norm(n5))
println("err=",norm(n6))
