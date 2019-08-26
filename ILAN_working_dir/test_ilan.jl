using NonlinearEigenproblems, Random, SparseArrays, Test, LinearAlgebra, PyPlot
import ..NEPSolver.ilan;
import ..NEPSolver.ilan_benchmark;
include("../src/method_ilan.jl");
include("../src/method_ilan_benchmark.jl");

n=20
LL=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)

x = range(0, stop = pi, length = n)
h=x[2]-x[1];
LL=LL/(h^2)
LL=-kron(LL,LL)
A=LL

b=broadcast((x,y)->sin(x+y),x,transpose(x))
B=sparse(1:n^2,1:n^2,b[:])

nep=DEP([A,B],[0,1.0])


v0=rand(n^2)

λ,v=ilan(nep,σ=0,γ=1;neigs=10,logger=1,maxit=100,tol=eps()*100,check_error_every=5,v=v0)
plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:black,linestyle=:none)
