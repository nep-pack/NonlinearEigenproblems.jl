using NonlinearEigenproblems, Random, SparseArrays, Test, LinearAlgebra, PyPlot, Revise, CSV
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


B=sparse(1:n^2,1:n^2,b[:])

nep=DEP([A,B],[0,1])
# define the relative error
nA1=opnorm(nep.A[1],Inf); nA2=opnorm(nep.A[2],Inf)
rel_err=(λ,v)->compute_resnorm(nep,λ,v)/((abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2)*norm(v));


v0=ones(n^2)

# COMPUTE REFERENCE EIGENVALUES WITH IAR
λ,v=tiar(nep;maxit=200,tol=1e-14,neigs=Inf,logger=1,check_error_every=Inf)
CSV.write("ILAN_figures/figure_1/dep_eigs.csv", λ)


# RUN ILAN (BEYN)
λ2,v2,err,_=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=50,tol=1e-10,check_error_every=1,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=-Inf,radius=4,N=1000))

CSV.write("ILAN_figures/figure_1/dep_ilan_Beyn.csv", λ2)

# EXPORT THE ERROR HIST
m,p=size(err);
for i=1:size(err,1) for j=1:size(err,2)	if err[i,j]==1 err[i,j]=NaN end end end
for j=1:p sort!(view(err,1:m,j);rev=true) end
for j=1:p semilogy(1:m,err[1:m,j],color="black",linestyle="-") end
writedlm( "ILAN_figures/figure_1/dep_ilan_Beyn_err.csv", [1:m err], ',')
#for j=1:p semilogy(1:m,err[1:m,j],color="black",linestyle="-") end


# RUN ILAN (PROJ)
λ3,v3,err3,_=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=50,tol=1e-10,check_error_every=1,v=v0,errmeasure=rel_err,
proj_solve=false)

CSV.write("ILAN_figures/figure_1/dep_ilan_Ritz.csv", λ3)

# EXPORT THE ERROR HIST
m,p=size(err3);
for i=1:size(err3,1) for j=1:size(err3,2)	if err3[i,j]==1 err3[i,j]=NaN end end end
for j=1:p sort!(view(err3,1:m,j);rev=true) end
for j=1:p semilogy(1:m,err3[1:m,j],color="black",linestyle="-") end
writedlm( "ILAN_figures/figure_1/dep_ilan_Ritz_err.csv", [1:m err3], ',')
#for j=1:p semilogy(1:m,err[1:m,j],color="black",linestyle="-") end
