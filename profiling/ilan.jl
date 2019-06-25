using NonlinearEigenproblems, Random, SparseArrays, Revise
include("../src/inner_solver.jl");
include("../src/method_ilan.jl");

n=4;

A1=[ 1 -4 -2
    -4  6  3
    -2  3  7]
A2=[ 3 5 4
     5 1 -1
     4 -1 11]
A3=[ 2   1 -11
     1   2  3
    -11  3 -3]
A4=[ -1  1   12
      1  -4  -2
     12  -2  7]

f1= S -> one(S)
f2= S -> -S
f3= S -> exp(-S)
f4= S -> sqrt(one(S)-2*S)


nep=SPMF_NEP([A1,A2,A3,A4],[f1,f2,f3,f4])
v0=ones(3)
V,H,ω=ilan(nep,σ=0,γ=1;neigs=10,v=v0,displaylevel=1,maxit=5,tol=eps()*100,check_error_every=1)
V
