using NonlinearEigenproblems, Random, SparseArrays, Revise, DelimitedFiles
import ..NEPSolver.ilan;
include("../src/method_ilan.jl");

# load the Voss symmetric DEP
n=200; nep=nep_gallery("dep_symm_double",n)

V,H,ω,HH=ilan(Float64,nep;Neig=10,displaylevel=1,maxit=200,tol=eps()*100,check_error_every=1)

# project (hardcoded for now)
Av=get_Av(nep); p = length(Av);
TT=promote_type(eltype(Av[1]),eltype(V))
pAv=Vector{Matrix{TT}}(undef, p)
for j=1:p
    pAv[j]=V'*(Av[j]*V)
end
pnep=DEP(pAv[2:p],nep.tauv)

err_lifted=(λ,z)->compute_resnorm(nep,λ,V*z)/n;

λ,_,err=tiar(pnep;Neig=200,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1,errmeasure=err_lifted)

m,p=size(err);

# sort error
for j=1:p err[1:m,j]=sort(err[1:m,j];rev=true) end

# plot conv hist
#for j=1:p semilogy(1:m,err[1:m,j],color="black",linestyle="-") end
#ylim(ymax=1)

# now export the error-matrix that will be loaded in tikz
err_print=ones(m,m+1)
err_print[1:m,1]=1:m
err_print[:,2:m+1]=err
writedlm("err_hist.csv",err_print,",")
