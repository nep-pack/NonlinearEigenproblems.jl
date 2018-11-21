using NonlinearEigenproblems, Random, SparseArrays, Revise, PyPlot
import ..NEPSolver.ilan;
include("../src/method_ilan.jl");

# n=1000
# Random.seed!(1) # reset the random seed
# K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
# A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
# A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
# A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';
# nep=DEP([A1,A2,A3],[0,1,0.8])

n=100
nep=nep_gallery("dep_symm_double",n)

σ=0;
γ=1;

#Dnep=DerSPMF(nep,σ,400)
V,H,ω,HH=ilan(Float64,nep;σ=σ,γ=γ,Neig=10,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)
#Q,_,_=svd(V)
Q=V;

# project (hardcoded for now)
Av=get_Av(nep); p = length(Av);
TT=promote_type(eltype(Av[1]),eltype(V))
pAv=Vector{Matrix{TT}}(undef, p)
for j=1:p
    pAv[j]=Q'*(Av[j]*Q)
end
pnep=DEP(pAv[2:p],nep.tauv)

#err_lifted=(λ,z)->compute_resnorm(nep,λ,Q*z)/n;
err_lifted=(λ,z)->compute_resnorm(nep,λ,Q*z)/n;

λ,_,err=iar(pnep;σ=σ,γ=γ,Neig=100,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1,errmeasure=err_lifted)

m,p=size(err);

# sort error
for j=1:p
    err[1:m,j]=sort(err[1:m,j];rev=true);
end

for j=1:p
    semilogy(1:m,err[1:m,j],color="black",linestyle="-");
end
ylim(ymax=1)
