using NonlinearEigenproblems, Random, SparseArrays, PyPlot, LinearAlgebra, CSV, DelimitedFiles


println("LOAD THE PROBLEM")
n=300; m=200
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
figure()
plot(real(λ1),imag(λ1),marker="*",markerfacecolor=:none,c=:black,linestyle=:none,label="eigenvalues (TIAR)")
