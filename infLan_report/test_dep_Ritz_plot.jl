using NonlinearEigenproblems, Random, SparseArrays, Revise, DelimitedFiles, PyPlot, LinearAlgebra
import ..NEPSolver.ilan;
import ..NEPSolver.iar;

include("../src/method_ilan.jl");
include("../src/method_iar.jl");
Random.seed!(1) # reset the random seed

# load the Voss symmetric DEP
n=1000
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, rand(3*n-2)); A1 = A1+A1';
A2 = sparse(K, J, rand(3*n-2)); A2 = A2+A2';
nep=DEP([A1,A2],[0,1])

# compute exact eigenvalues
λ1,_,err,V1,_,H1=iar(nep;Neig=100,displaylevel=1,maxit=150,tol=eps()*100,check_error_every=1)

m=100
V2,H2,ω2,_,_,WW=ilan(nep;Neig=10,displaylevel=1,maxit=m,tol=eps()*100,check_error_every=1)
λ2=1 ./ eigvals(H2[1:m,1:m])
#V,_,_=svd(V2)    # manually reorth
V=V2
V=V[:,1:m]

# Create a projected NEP
mm=size(V,2)
pnep=create_proj_NEP(nep,mm); # maxsize=mm
set_projectmatrices!(pnep,V,V);
λ3,_,err,V3,_,H3=iar(pnep;Neig=200,displaylevel=1,maxit=100,tol=eps()*1000,check_error_every=1)

# plot the eigenvalues and Ritz values
figure()
scatter(real(λ1),imag(λ1),marker="o",c=:red)        # original eigenvalues
scatter(real(λ2),imag(λ2),marker="d",c=:none)       # Ritz values
scatter(real(λ3),imag(λ3),marker="s",c=:none)       # eigenvalues of proj problem

# define the error function

nA1=opnorm(nep.A[1],Inf); nA2=opnorm(nep.A[2],Inf)
rel_err=(λ,z)->compute_resnorm(nep,λ,z)/(abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2);

# compute and plot the different convergence histories
# for ilan
err1=ones(m,m)
for j=1:m
    F=eigen(H2[1:j,1:j])
    λ=1 ./ F.values; z=F.vectors;
    for i=1:j
        err1[j,i]=rel_err(λ[i],WW[:,1:j]*z[:,i]);
    end
end

# for the projected problem
err2=ones(m,m)
for j=1:m
    F=eigen(H3[1:j,1:j])
    λ=1 ./ F.values; z=F.vectors;
    for i=1:j
        err2[j,i]=rel_err(λ[i],V*(V3[1:m,1:j]*z[:,i]));
    end
end

figure()
# sort error
for j=1:m err1[1:m,j]=sort(err1[1:m,j];rev=true) end
for j=1:m err2[1:m,j]=sort(err2[1:m,j];rev=true) end

# plot conv hist
for j=1:m semilogy(1:m,err1[1:m,j],color="red",linestyle="-") end
for j=1:m semilogy(1:m,err2[1:m,j],color="black",linestyle="-") end
ylim(ymax=1)
1
