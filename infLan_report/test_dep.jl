using NonlinearEigenproblems, Random, SparseArrays, Revise, PyPlot
import ..NEPSolver.ilan;
include("../src/method_ilan.jl");

println("Loading the DEP:", `dep_symm_double`)
n=20; nep=nep_gallery("dep_symm_double",n)
# define the relative error
nA1=opnorm(nep.A[1],Inf); nA2=opnorm(nep.A[2],Inf)
rel_err=(λ,v)->compute_resnorm(nep,λ,v)/((abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2)*norm(v));

println("Run the infinite Lanczos method")
λ,W,V=ilan(Float64,nep;Neig=10,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=Inf,errmeasure=rel_err)

println("Number of computed eigenpairs: ", length(λ))
for j=1:length(λ)
    println("Residual of the eigepair ", j, "th = ",rel_err(λ[j],W[:,j]))
end

# Create a projected NEP
mm=size(V,2)
pnep=create_proj_NEP(nep,mm);
set_projectmatrices!(pnep,V,V);

err_lifted=(λ,z)->rel_err(λ,V*z);
println("At the last iteration re-solve the projected problem with", `iar`,"to see the convergence history")
λ,_,err=iar(pnep;Neig=200,displaylevel=0,maxit=100,tol=eps()*100,check_error_every=1,errmeasure=err_lifted)

m,p=size(err);
for j=1:p sort!(view(err,1:m,j);rev=true) end
for j=1:p semilogy(1:m,err[1:m,j],color="black",linestyle="-") end
ylim(ymax=10)
