using NonlinearEigenproblems, Random, SparseArrays, Revise, PyPlot
import ..NEPSolver.ilan;
include("../src/method_ilan.jl");

# load the Voss symmetric DEP
n=100; nep=nep_gallery("dep_symm_double",n)

V,H,ω,HH=ilan(Float64,nep;Neig=10,displaylevel=1,maxit=200,tol=eps()*100,check_error_every=1)

# Create a projected NEP
mm=size(V,2)
pnep=create_proj_NEP(nep,mm); # maxsize=mm
set_projectmatrices!(pnep,V,V);

err_lifted=(λ,z)->compute_resnorm(nep,λ,V*z)/n;

λ,_,err=tiar(pnep;Neig=200,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1,errmeasure=err_lifted)

m,p=size(err);

# sort error
for j=1:p err[1:m,j]=sort(err[1:m,j];rev=true) end

# plot conv hist
for j=1:p semilogy(1:m,err[1:m,j],color="black",linestyle="-") end
ylim(ymax=1)
