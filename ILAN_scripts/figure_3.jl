using NonlinearEigenproblems, Random, SparseArrays, PyPlot, LinearAlgebra

pygui(true)


println("LOAD THE PROBLEM")
n=500; m=200
A1=spdiagm(-1 => -ones(n-1), 0 => -2*ones(n), 1 => -ones(n-1))*n
A2=one(A1)
A3=spdiagm(-1 => ones(n-1), 0 => 0*ones(n) , 1 => ones(n-1))
A4=spdiagm(-1 =>  1im*ones(n-1) )
f1= S -> one(S); f2= S -> -S; f3= S -> sin(S); f4= S -> exp(-S)
nep1=SPMF_NEP([A1,A2,A3,A4],[f1,f2,f3,f4])

nA1=norm(A1); nA2=norm(A2); nA3=norm(A3); nA4=norm(A4);
rel_err=(λ,z)->compute_resnorm(nep1,λ,z)/((nA1*abs(f1(λ))+nA2*abs(f2(λ))+nA3*abs(f3(λ))+nA4*abs(f4(λ)))*norm(z));

println("COMPUTE EIGENVALUES WITH TIAR")
λ1,_,err_iar=tiar(nep1;neigs=Inf,logger=0,maxit=m,tol=1e-12,check_error_every=1,errmeasure=rel_err)


println("CONSTRUCTING THE SYMMETRIZED PROBLEM")
O=zero(A1);
AA1=[O A1; transpose(A1) O];
AA2=[O A2; transpose(A2) O];
AA3=[O A3; transpose(A3) O];
AA4=[O A4; transpose(A4) O];

nep=SPMF_NEP([AA1,AA2,AA3,AA4],[f1,f2,f3,f4])
nep=DerSPMF(nep,0,2*m)

println("SOLVING THE PROBLEM WITH ILAN")
Σ=float([0.5+0.5im, -1.5-0.5im, -1.5-1.5im,0.5-1.5im])

v0=ones(size(nep,1))
λ,W,err,_=ilan(nep;v=v0,neigs=Inf,logger=1,maxit=50,tol=1e-8,check_error_every=1,errmeasure=(λ,v)->rel_err(λ,v[n+1:end]),inner_solver_method=NEPSolver.NleigsInnerSolver(Σ=Σ,tol=1e-2))
W = W[n+1:end,:]    # extract the eigenvectors

# EXPORT THE ERROR HIST
m,p=size(err);
for i=1:size(err,1) for j=1:size(err,2)	if err[i,j]==1 err[i,j]=NaN end end end
for j=1:p sort!(view(err,1:m,j);rev=true) end

figure(1)
m,p=size(err);
for j=1:p semilogy(1:m,err[1:m,j],color="black",linestyle="-") end
ylim(ymax=10)

figure(2)
plot(real(λ1),imag(λ1),marker="*",markerfacecolor=:none,c=:black,linestyle=:none,label="eigenvalues (TIAR)")
plot(real(λ),imag(λ),marker="o",markerfacecolor=:none,c=:red,linestyle=:none,label="ILAN.")
plot([-1.5,-1.5,0.5,0.5,-1.5],[-1.5,0.5,0.5,-1.5,-1.5],linestyle="-",color="black", linewidth=2.0,label="region of interest")
legend()
legend(numpoints = 1)	# display only one marker
