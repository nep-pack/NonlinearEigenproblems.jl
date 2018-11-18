using NonlinearEigenproblems, Random, SparseArrays, Revise, PyPlot
import ..NEPSolver.ilan;
include("../src/method_ilan.jl");

n=1000
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
A3 = sparse(K, J, rand(3*n-2)); A3 = A3+A3';

nep=DEP([A1,A2,A3],[0,1,0.8])

σ=0;
γ=1;

Dnep=DerSPMF(nep,σ,400)
V,H,ω,HH=ilan(Dnep;σ=σ,γ=γ,Neig=10,displaylevel=1,maxit=200,tol=eps()*100,check_error_every=1)
#Q,_,_=svd(V)
Q=V;

# project (hardcoded for now)
AA1=Q'*(A1*Q);
AA2=Q'*(A2*Q);
AA3=Q'*(A3*Q);

#err_lifted=(λ,z)->compute_resnorm(nep,λ,Q*z)/n;
err_lifted=(λ,z)->compute_resnorm(nep,λ,Q*z)/n;

pnep=DEP([AA1,AA2,AA3],[0,1,0.8])
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
