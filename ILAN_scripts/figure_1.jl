using NonlinearEigenproblems, Random, SparseArrays, LinearAlgebra, PyPlot


n=20 # FIX TO 100

pygui(true)

# DISCRETIZE LAPLACIAN
LL=-sparse(1:n,1:n,2*ones(n))+sparse(2:n,1:n-1,ones(n-1),n,n)+sparse(1:n-1,2:n,ones(n-1),n,n)
x = range(0, stop = π, length = n);
h=x[2]-x[1]; LL=LL/(h^2); A=-kron(LL,LL)
# DISCRETIZE OTHER LINEAR TERM
b=broadcast((x,y)->-x*sin(y-x),x,transpose(x))
B=sparse(1:n^2,1:n^2,b[:])
# CONSTRUCT THE DEP
nep=DEP([A,B],[0,1])
# DEFINE THE RELAIVE ERROR
nA1=opnorm(nep.A[1],Inf); nA2=opnorm(nep.A[2],Inf)
rel_err=(λ,v)->compute_resnorm(nep,λ,v)/((abs(λ)+abs(exp(nep.tauv[1]))*nA1+abs(exp(nep.tauv[2]))*nA2)*norm(v));


v0=ones(n^2)

# COMPUTE REFERENCE EIGENVALUES WITH IAR
λ,v=tiar(nep;maxit=200,tol=1e-14,neigs=Inf,logger=1,check_error_every=Inf)

# RUN ILAN (BEYN)
λ2,v2,err2,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=1,maxit=50,tol=1e-8,check_error_every=1,v=v0,errmeasure=rel_err,
inner_solver_method=NEPSolver.ContourBeynInnerSolver(tol=1e-8,radius=4,N=1000))

# PLOT ERROR HIST
m,p=size(err2);
for i=1:size(err2,1) for j=1:size(err2,2)	if err2[i,j]==1 err2[i,j]=NaN end end end
for j=1:p sort!(view(err2,1:m,j);rev=true) end
for j=1:p semilogy(1:m,err2[1:m,j],color="black",linestyle="-") end

# RUN ILAN (PROJ)
λ3,v3,err3,_=ilan(nep,σ=0,γ=1;neigs=Inf,logger=1,maxit=50,tol=1e-8,check_error_every=1,v=v0,errmeasure=rel_err,
proj_solve=false)

# PLOT ERROR HIST
m,p=size(err3);
for i=1:size(err3,1) for j=1:size(err3,2)	if err3[i,j]==1 err3[i,j]=NaN end end end
for j=1:p sort!(view(err3,1:m,j);rev=true) end
for j=1:p semilogy(1:m,err3[1:m,j],color="red",linestyle="--") end


# PLOT THE SPECTRUM AND THE CONVERGED APPROXIMATIONS
#plot(real(λ2),imag(λ2),marker="o",markerfacecolor=:none,c=:red,linestyle=:none)
