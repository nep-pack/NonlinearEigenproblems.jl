using NonlinearEigenproblems, Random, SparseArrays, Test, LinearAlgebra, PyPlot, Revise, CSV, DelimitedFiles

import ..NEPSolver.ilan;
include("../src/method_ilan.jl");


# construct the Laplacian
n=100
LL=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)

x = range(0, stop = π, length = n)
h=x[2]-x[1];
LL=LL/(h^2)
LL=-kron(LL,LL)

# construct the other matrix coefficient
b=broadcast((x,y)->-x*sin(y-x),x,transpose(x))

A=LL; B=sparse(1:n^2,1:n^2,b[:])

nep=DEP([A,B],[0,1.0])

# define the relative error
nA1=opnorm(nep.A[1],Inf); nA2=opnorm(nep.A[2],Inf)
rel_err=(λ,v)->compute_resnorm(nep,λ,v)/((abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2)*norm(v));

v0=ones(n^2)
tol_used=1e-6
# COMPUTE REFERENCE EIGENVALUES WITH IAR
println("COMPUTING REFERENCE EIGENVALUES")
λ,v=tiar(nep;maxit=350,tol=1e-10,neigs=Inf,logger=1,check_error_every=Inf,errmeasure=rel_err
)
CSV.write("ILAN_figures/figure_1/dep_eigs.csv", λ)

# COMPUTE REFERENCE EIGENVALUES WITH ILAN (Beyn)
println("RUNNING ILAN (BEYN)")
λ2,v2,err2,_=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=50,tol=tol_used,check_error_every=1,v=v0,
inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=-Inf,radius=4,N=1000),errmeasure=rel_err
)
# save computed eigenvalues
CSV.write("ILAN_figures/figure_1/dep_ilan_Beyn.csv", λ2)
# save err hist
m,p=size(err2);
for j=1:p sort!(view(err2,1:m,j);rev=true) end
for i=1:size(err2,1) for j=1:size(err2,2) if err2[i,j]==1 err2[i,j]=NaN end end end
writedlm( "ILAN_figures/figure_1/dep_ilan_Beyn_err.csv", err2, ',')

# COMPUTE REFERENCE EIGENVALUES EXTRACTING RITZ VALUES
println("RUNNING ILAN (RITZ)")
λ3,v3,err3,_=ilan(nep,σ=0,γ=1;neigs=100,logger=1,maxit=50,tol=tol_used,check_error_every=1,v=v0,
proj_solve=false,errmeasure=rel_err)
# save computed eigenvalues
CSV.write("ILAN_figures/figure_1/dep_ilan_Ritz.csv", λ3)
# save err hist
m,p=size(err3);
for j=1:p sort!(view(err3,1:m,j);rev=true) end
for i=1:size(err3,1) for j=1:size(err3,2) if err3[i,j]==1 err3[i,j]=NaN end end end
writedlm( "ILAN_figures/figure_1/dep_ilan_Ritz_err.csv", err3, ',')
