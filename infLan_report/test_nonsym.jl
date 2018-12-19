using NonlinearEigenproblems, Random, SparseArrays, Revise, PyPlot, DelimitedFiles
import ..NEPSolver.ilan;
import ..NEPSolver.iar;

include("../src/method_ilan.jl");
include("../src/method_iar.jl");

n=500; m=250

A1=spdiagm(-1 => -ones(n-1), 0 => -2*ones(n), 1 => -ones(n-1))*n
A2=one(A1)
A3=spdiagm(-1 => ones(n-1), 0 => 0*ones(n) , 1 => ones(n-1))
A4=spdiagm(-1 =>  1im*ones(n-1) )

f1= S -> one(S)
f2= S -> -S
f3= S -> S*sin(S)
f4= S -> exp(-S)
nep1=SPMF_NEP([A1,A2,A3,A4],[f1,f2,f3,f4])

nA1=norm(A1); nA2=norm(A2); nA3=norm(A3); nA4=norm(A4);

rel_err=(λ,z)->compute_resnorm(nep1,λ,z)/((
nA1*abs(f1(λ))+
nA2*abs(f2(λ))+
nA3*abs(f3(λ))+
nA4*abs(f4(λ))
)*norm(z));

λ1,_,err_iar=tiar(nep1;Neig=100,displaylevel=1,maxit=m,tol=1e-12,check_error_every=1,errmeasure=rel_err)

O=zero(A1);
AA1=[O A1; transpose(A1) O];
AA2=[O A2; transpose(A2) O];
AA3=[O A3; transpose(A3) O];
AA4=[O A4; transpose(A4) O];

nep=SPMF_NEP([AA1,AA2,AA3,AA4],[f1,f2,f3,f4])
nep=DerSPMF(nep,0,2*m)

λ,W,_=ilan(nep;Neig=200,displaylevel=1,maxit=100,tol=1e-6,check_error_every=Inf,errmeasure=(λ,v)->rel_err(λ,v[n+1:end]))
W = W[n+1:end,:]    # extract the eigenvectors

println("Number of computed eigenpairs: ", length(λ))
for j=1:length(λ)
    println("Residual of the eigepair ", j, "th = ",rel_err(λ[j],W[:,j]))
end

plot(real(λ1),imag(λ1),marker="o",markerfacecolor=:none,c=:red,linestyle=:none,label="eigenvalues (TIAR)")
plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:black,linestyle=:none,label="INF. LAN.")
legend()
