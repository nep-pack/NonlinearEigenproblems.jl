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
tol_used=1e-10
# COMPUTE REFERENCE EIGENVALUES WITH IAR
println("COMPUTING REFERENCE EIGENVALUES")
λ,v=tiar(nep;maxit=350,tol=1e-10,neigs=Inf,logger=1,check_error_every=Inf,errmeasure=rel_err
)
CSV.write("ILAN_figures/dep_eigs.csv", λ)

# COMPUTE REFERENCE EIGENVALUES WITH ILAN (Beyn)
println("RUNNING ILAN (BEYN)")
λ2,v2=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=100,tol=tol_used,check_error_every=Inf,v=v0,
inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=1e-2,radius=3,N=1000),errmeasure=rel_err
)
CSV.write("ILAN_figures/dep_ilan_Beyn.csv", λ2)

# COMPUTE REFERENCE EIGENVALUES WITH ILAN (nleigs)
println("RUNNING ILAN (NLEIGS)")
θ=range(0,stop=2π,length=1000); r=6; Σ=r*cos.(θ) + 1im*r*sin.(θ)
λ3,v3=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=100,tol=tol_used,check_error_every=Inf,v=v0,
inner_solver_method=NEPSolver.NleigsInnerSolver(Σ=Σ,tol=1e-2),errmeasure=rel_err
)
CSV.write("ILAN_figures/dep_ilan_nleigs.csv", λ3)

# COMPUTE REFERENCE EIGENVALUES WITH ILAN (iar)
println("RUNNING ILAN (IAR)")
λ4,v4=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=100,tol=tol_used,check_error_every=Inf,v=v0,
inner_solver_method=NEPSolver.IARInnerSolver(tol=1e-2),errmeasure=rel_err
)
CSV.write("ILAN_figures/dep_ilan_iar.csv", λ4)

# COMPUTE REFERENCE EIGENVALUES EXTRACTING RITZ VALUES
println("RUNNING ILAN (RITZ)")
λ5,v5=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=100,tol=Inf,check_error_every=Inf,v=v0,
proj_solve=false,errmeasure=rel_err)

CSV.write("ILAN_figures/dep_ilan_Ritz.csv", λ5)
