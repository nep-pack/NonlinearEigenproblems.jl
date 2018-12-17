using NonlinearEigenproblems, Random, SparseArrays, Revise, PyPlot
import ..NEPSolver.ilan;
import ..NEPSolver.iar;

include("../src/method_ilan.jl");
include("../src/method_iar.jl");
mm=10

nep1=nep_gallery("schrodinger_movebc")
Av=get_Av(nep1); fv=get_fv(nep1);
A1=Av[1]+0im*Av[1];   f1=fv[1];
A2=Av[2]+0im*Av[1];   f2=fv[2];
A3=Av[3]+0im*Av[1];   f3=fv[3];
A4=Av[4]+0im*Av[1];   f4=fv[4];
nep1=SPMF_NEP([A1,A2,A3,A4],[f1,f2,f3,f4])
nep1=DerSPMF(nep1,σ,mm)



nA1=norm(A1); nA2=norm(A2); nA3=norm(A3); nA4=norm(A4);

rel_err=(λ,z)->compute_resnorm(nep1,λ,z)/(
nA1*abs(f1(λ))+
nA2*abs(f2(λ))+
nA3*abs(f3(λ))+
nA4*abs(f4(λ))
);
σ=0;
γ=1;
λ,_,err=iar(nep1;σ=σ,γ=γ,Neig=100,displaylevel=1,maxit=mm,tol=1e-12,check_error_every=1,errmeasure=rel_err)


O=zero(A1);
AA1=[O A1; transpose(A1) O];
AA2=[O A2; transpose(A2) O];
AA3=[O A3; transpose(A3) O];
AA4=[O A4; transpose(A4) O];


nep=SPMF_NEP([AA1,AA2,AA3,AA4],[f1,f2,f3,f4])
nep=DerSPMF(nep,σ,mm)

V,H,ω,HH=ilan(nep;σ=σ,γ=γ,Neig=200,displaylevel=1,maxit=mm,tol=eps()*100,check_error_every=1)
V,_,_=svd(V)
Q=V;

# Create a projected NEP
mm=size(V,2)
pnep=create_proj_NEP(nep,mm); # maxsize=mm
set_projectmatrices!(pnep,V,V);


#err_lifted=(λ,z)->compute_resnorm(nep,λ,Q*z)/n;
err_lifted=(λ,z)->compute_resnorm(nep1,λ,(Q*z)[n+1:end])/(
nA1*abs(f1(λ))+
nA2*abs(f2(λ))+
nA3*abs(f3(λ))+
nA4*abs(f4(λ))
);
λ1,_,err=iar(pnep;σ=σ,γ=γ,Neig=100,displaylevel=1,maxit=100,tol=1e-9,check_error_every=1,errmeasure=err_lifted)

m,p=size(err);

# sort error
for j=1:p
    err[1:m,j]=sort(err[1:m,j];rev=true);
end

for j=1:p
    semilogy(1:m,err[1:m,j],color="black",linestyle="-");
end
ylim(ymax=1)

figure()
plot(real(λ1),imag(λ1),marker="o",markerfacecolor=:none,c=:red,linestyle=:none)         # Ritz values
plot(real(λ),imag(λ),marker="*",markerfacecolor=:none,c=:black,linestyle=:none)         # Ritz values

# plot err hist in lest and right eigenvectors
