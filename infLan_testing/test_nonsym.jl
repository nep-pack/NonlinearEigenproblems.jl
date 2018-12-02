using NonlinearEigenproblems, Random, SparseArrays, Revise, PyPlot
import ..NEPSolver.ilan;
import ..NEPSolver.iar;

include("../src/method_ilan.jl");
include("../src/method_iar.jl");


n=300
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2));
A2 = sparse(K, J, rand(3*n-2));
A3 = sparse(K, J, rand(3*n-2));
A4 = sparse(K, J, rand(3*n-2));

f1= S -> one(S)
f2= S -> -S
f3= S -> exp(-S)
nep1=SPMF_NEP([A1,A2,A3],[f1,f2,f3])

A1=Matrix(A1); A2=Matrix(A2);
A3=Matrix(A3); A4=Matrix(A4);

O=zero(A1);
A1=[O A1; A1' O];
A2=[O A2; A2' O];
A3=[O A3; A3' O];
A4=[O A4; A4' O];

nep=SPMF_NEP([A1,A2,A3],[f1,f2,f3])

σ=0;
γ=1;

#Dnep=DerSPMF(nep,σ,400)
V,H,ω,HH=ilan(Float64,nep;σ=σ,γ=γ,Neig=200,displaylevel=1,maxit=200,tol=eps()*100,check_error_every=1)
V,_,_=svd(V)
#Q=V;

# Create a projected NEP
mm=size(V,2)
pnep=create_proj_NEP(nep,mm); # maxsize=mm
set_projectmatrices!(pnep,V,V);


#err_lifted=(λ,z)->compute_resnorm(nep,λ,Q*z)/n;
err_lifted=(λ,z)->compute_resnorm(nep1,λ,(Q*z)[n+1:end])/n;
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
