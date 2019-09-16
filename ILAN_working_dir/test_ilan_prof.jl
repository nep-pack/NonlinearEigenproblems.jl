using NonlinearEigenproblems, Random, SparseArrays, Test, LinearAlgebra, PyPlot, Profile, BenchmarkTools
import ..NEPSolver.ilan;
import ..NEPSolver.iar;

#Profile.init()
#Profile.init(n = 10^7, delay = 0.1)
include("../src/method_ilan.jl");
include("../src/method_iar.jl");


n=100
LL=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)

x = range(0, stop = π, length = n)
h=x[2]-x[1];
LL=LL/(h^2)
LL=-kron(LL,LL)
A=LL

b=broadcast((x,y)->-x*sin(y-x),x,transpose(x))
#b=broadcast((x,y)->-sin(y-x),x,transpose(x))
#b=broadcast((x,y)->sin(x^(π-y)),x,transpose(x))

B=sparse(1:n^2,1:n^2,b[:])

nep=DEP([A,B],[0,1])
# define the relative error
nA1=opnorm(nep.A[1],Inf); nA2=opnorm(nep.A[2],Inf)
rel_err=(λ,v)->compute_resnorm(nep,λ,v)/((abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2)*norm(v));


v0=ones(n^2)

@btime λ2,v2,err,_=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=50,tol=1e-10,check_error_every=Inf,v=v0,errmeasure=rel_err,inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=-Inf,radius=4,N=1000))
#λ2,v2,err,_=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=50,tol=1e-10,check_error_every=Inf,v=v0,errmeasure=rel_err,inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=-Inf,radius=4,N=1000))
#@profile λ2,v2,err,_=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=50,tol=1e-10,check_error_every=Inf,v=v0,errmeasure=rel_err,inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=-Inf,radius=4,N=1000))

#Profile.print(format=:flat)
1
