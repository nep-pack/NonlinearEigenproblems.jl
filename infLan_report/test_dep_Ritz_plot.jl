using NonlinearEigenproblems, Random, SparseArrays, Revise, DelimitedFiles, PyPlot, LinearAlgebra
import ..NEPSolver.ilan;
import ..NEPSolver.iar;

include("../src/method_ilan.jl");
include("../src/method_iar.jl");
Random.seed!(1) # reset the random seed

# load the Voss symmetric DEP
n=300
Random.seed!(1) # reset the random seed
K = [1:n;2:n;1:n-1]; J=[1:n;1:n-1;2:n]
A1 = sparse(K, J, n*rand(3*n-2)); A1 = A1'+A1;
A2 = sparse(K, J, n*rand(3*n-2)); A2 = A2'+A2;

nep=DEP([A1,A2],[0,1])

# define the error function
nA1=opnorm(nep.A[1],Inf); nA2=opnorm(nep.A[2],Inf)
rel_err=(λ,z)->compute_resnorm(nep,λ,z)/(abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2);

# compute exact eigenvalues
λ1,_,err,V1,_,H1=iar(nep;Neig=100,displaylevel=1,maxit=200,tol=eps()*100,check_error_every=1,errmeasure=rel_err)

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
mm=100;

# filter converged Ritz pairs
F=eigen(H2[1:mm,1:mm])
λt=1 ./ F.values; z=F.vectors;
λ2=[]
for i=1:mm
    if rel_err(λt[i],WW[:,1:mm]*z[:,i])<1e-12
        push!(λ2,λt[i])
    end
end

F=eigen(H3[1:mm,1:mm])
λt=1 ./ F.values; z=F.vectors;
λ3=[]
for i=1:mm
    if rel_err(λt[i],V*(V3[1:mm,1:mm]*z[:,i]))<1e-12
        push!(λ3,λt[i])
    end
end

plot(real(λ1),imag(λ1),marker="+",markerfacecolor=:none,c=:black,linestyle=:none)       # original eigenvalues
plot(real(λ2),imag(λ2),marker="d",markerfacecolor=:none,c=:red,linestyle=:none)         # Ritz values
plot(real(λ3),imag(λ3),marker="o",markerfacecolor=:none,c=:blue,linestyle=:none)        # eigenvalues of proj NonlinearEigenproblems

# save in a file the eigenvalues and Ritz pairs
writedlm("Ritz_comp_l1.csv",[real(λ1) imag(λ1)],",")
writedlm("Ritz_comp_l2.csv",[real(λ2) imag(λ2)],",")
writedlm("Ritz_comp_l3.csv",[real(λ3) imag(λ3)],",")

# compute and plot the different convergence histories
# for ilan
err1=NaN*ones(m,m)
for j=1:m
    F=eigen(H2[1:j,1:j])
    λ=1 ./ F.values; z=F.vectors;
    for i=1:j
        err1[j,i]=rel_err(λ[i],WW[:,1:j]*z[:,i]);
    end
end

# for the projected problem
err2=NaN*ones(m,m)
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

# now export the error-matrix that will be loaded in tikz
err_print=NaN*ones(m,m+1)
err_print[1:m,1]=1:m
err_print[:,2:m+1]=err1
writedlm("Ritz_comp_err1.csv",err_print,",")

err_print=NaN*ones(m,m+1)
err_print[1:m,1]=1:m
err_print[:,2:m+1]=err2
writedlm("Ritz_comp_err2.csv",err_print,",")
