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
