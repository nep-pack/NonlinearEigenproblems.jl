using NonlinearEigenproblems, Random, SparseArrays, Revise, DelimitedFiles, PyPlot
import ..NEPSolver.ilan;
import ..NEPSolver.tiar;
include("../src/method_ilan.jl");
include("../src/method_tiar.jl");

# load the Voss symmetric DEP
n=320; nep=nep_gallery("dep_symm_double",n)

V,H,ω,HH=ilan(Float64,nep;Neig=10,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)
V,_,_=svd(V)    # manually reorth

# Create a projected NEP
mm=size(V,2)
pnep=create_proj_NEP(nep,mm); # maxsize=mm
set_projectmatrices!(pnep,V,V);

nA1=opnorm(nep.A[1],Inf)
nA2=opnorm(nep.A[2],Inf)
err_lifted=(λ,z)->compute_resnorm(nep,λ,V*z)/(abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2);
λ,_,err=tiar(pnep;Neig=200,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1,errmeasure=err_lifted)

m,p=size(err);

# sort error
for j=1:p err[1:m,j]=sort(err[1:m,j];rev=true) end

# plot conv hist
for j=1:p semilogy(1:m,err[1:m,j],color="black",linestyle="-") end
ylim(ymax=1)

# now export the error-matrix that will be loaded in tikz
err_print=NaN*ones(m,m+1)
err_print[1:m,1]=1:m
err_print[:,2:m+1]=err
writedlm("err_hist.csv",err_print,",")
