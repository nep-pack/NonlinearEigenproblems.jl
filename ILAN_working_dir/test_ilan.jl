using NonlinearEigenproblems, Random, SparseArrays, Test, LinearAlgebra, PyPlot, Revise
import ..NEPSolver.ilan;
import ..NEPSolver.ilan_benchmark;

include("../src/method_ilan.jl");
include("../src/method_ilan_benchmark.jl");


n=50
LL=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)

x = range(0, stop = pi, length = n)
h=x[2]-x[1];
LL=LL/(h^2)
LL=-kron(LL,LL)
A=LL

b=broadcast((x,y)->sin(x+y)+1,x,transpose(x))
B=sparse(1:n^2,1:n^2,b[:])

nep=DEP([A,B],[0,1.0])
# define the relative error
nA1=opnorm(nep.A[1],Inf); nA2=opnorm(nep.A[2],Inf)
rel_err=(λ,v)->compute_resnorm(nep,λ,v)/((abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2)*norm(v));


v0=rand(n^2)


# COMPUTE REFERENCE EIGENVALUES WITH IAR
#λ,v=iar(nep;maxit=100,tol=1e-6,neigs=Inf,logger=1)
plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:black,linestyle=:none)

# COMPUTE EIGENVALUES WITH
λ2,v2=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=100,tol=1e-6,check_error_every=20,v=v0,inner_solver_method=IARInnerSolver)
plot(real(λ2),imag(λ2),marker="o",markerfacecolor=:none,c=:black,linestyle=:none)

# todo: fix too many unconv eigs displayed
