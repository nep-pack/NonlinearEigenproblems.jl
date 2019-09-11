using NonlinearEigenproblems, Random, SparseArrays, Test, LinearAlgebra, PyPlot, Revise
import ..NEPSolver.ilan;
import ..NEPSolver.ilan_benchmark;
import ..NEPSolver.iar;


include("../src/method_ilan.jl");
include("../src/method_ilan_benchmark.jl");
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


# COMPUTE REFERENCE EIGENVALUES WITH IAR
λ,v=tiar(nep;maxit=100,tol=1e-14,neigs=Inf,logger=1,check_error_every=Inf)
#plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:black,linestyle=:none)
# θ=range(0,stop=2π,length=1000); r=6; Σ=r*cos.(θ) + 1im*r*sin.(θ)
#
#COMPUTE EIGENVALUES WITH ILAN
λ2,v2,err,_=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=50,tol=1e-10,check_error_every=1,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=-Inf,radius=4,N=1000))
#inner_solver_method=NEPSolver.IARInnerSolver(tol=1e2,maxit=50))
#inner_solver_method=NEPSolver.NleigsInnerSolver(Σ=Σ,tol=1e-2))

#λ2,v2=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=100,tol=1e-6,check_error_every=1,v=v0,proj_solve=false)

#plot(real(λ2),imag(λ2),marker="o",markerfacecolor=:none,c=:black,linestyle=:none)

err1=err
pygui(true)
m,p=size(err);
for j=1:p sort!(view(err,1:m,j);rev=true) end
for j=1:p semilogy(1:m,err[1:m,j],color="black",linestyle="-") end
ylim(ymax=10)

figure()
plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:black,linestyle=:none)
plot(real(λ2),imag(λ2),marker="o",markerfacecolor=:none,c=:black,linestyle=:none)
