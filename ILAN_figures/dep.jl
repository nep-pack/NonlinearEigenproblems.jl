using NonlinearEigenproblems, Random, SparseArrays, Test, LinearAlgebra, PyPlot, Revise, CSV
import ..NEPSolver.ilan;
import ..NEPSolver.ilan_benchmark;
import ..NEPSolver.iar;


include("../src/method_ilan.jl");
include("../src/method_ilan_benchmark.jl");
include("../src/method_iar.jl");


n=100
LL=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)

x = range(0, stop = pi, length = n)
h=x[2]-x[1];
LL=LL/(h^2)
LL=-kron(LL,LL)
A=LL

b=broadcast((x,y)->sin(x+y),x,transpose(x))
B=sparse(1:n^2,1:n^2,b[:])

nep=DEP([A,B],[0,1.0])
# define the relative error
nA1=opnorm(nep.A[1],Inf); nA2=opnorm(nep.A[2],Inf)
rel_err=(λ,v)->compute_resnorm(nep,λ,v)/((abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2)*norm(v));


v0=rand(n^2)

# COMPUTE REFERENCE EIGENVALUES WITH IAR
λ,v=tiar(nep;maxit=250,tol=1e-12,neigs=Inf,logger=1,check_error_every=Inf)
plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:black,linestyle=:none)
CSV.write("ILAN_figures/dep_eigs.csv", λ)
