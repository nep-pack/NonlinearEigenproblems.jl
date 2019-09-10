using NonlinearEigenproblems, Random, SparseArrays, Test, LinearAlgebra, PyPlot, Revise
import ..NEPSolver.ilan;
import ..NEPSolver.ilan_benchmark;
import ..NEPSolver.iar;


include("../src/method_ilan.jl");
include("../src/method_ilan_benchmark.jl");
include("../src/method_iar.jl");


n=50
LL=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)

x = range(0, stop = pi, length = n)
h=x[2]-x[1];



A0=sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)
A0=A0/h^2

A1=sparse(1:n^2,1:n^2,a0[:])
A2=sparse(1:n^2,1:n^2,a1[:])


nep=DEP([A,B],[0,1.0])
# define the relative error
nA1=opnorm(nep.A[1],Inf); nA2=opnorm(nep.A[2],Inf)
rel_err=(λ,v)->compute_resnorm(nep,λ,v)/((abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2)*norm(v));


v0=ones(n^2)

# COMPUTE REFERENCE EIGENVALUES WITH IAR
@time λ1,v1=iar(nep;maxit=100,tol=1e-10,neigs=Inf,logger=0,check_error_every=Inf)
#@time λ2,v2=tiar(nep;maxit=100,tol=1e-6,neigs=Inf,logger=0,check_error_every=Inf)
@time λ3,v3=ilan(nep,σ=0,γ=1;neigs=100,logger=0,maxit=100,tol=1e-10,check_error_every=Inf,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=200))
1
