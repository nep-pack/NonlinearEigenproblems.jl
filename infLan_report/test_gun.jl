using NonlinearEigenproblems, Random, SparseArrays, Revise, PyPlot
import ..NEPSolver.ilan;
include("../src/method_ilan.jl");
nep = nep_gallery("nlevp_native_gun")

# manually loading the matrices
Av=get_Av(nep)
K=Av[1]; M=-Av[2]; W1=Av[3]; W2=Av[4]
nK=opnorm(K,1); nM=opnorm(M,1); nW1=opnorm(W1,1); nW2=opnorm(W2,1);


# define the functions
σ=108.8774; α=(300^2-200^2)/10;  λ0=250^2; # scale and shift
a1 =α; b1=λ0; a2=α; b2=λ0-σ^2
f1 = l-> one(l);
f2 = l-> l;
f3 = l-> sqrt(a1*l+b1*one(l));
f4 = l-> sqrt(a2*l+b2*one(l));

nep=SPMF_NEP([K-λ0*M,-α*M,1im*W1,1im*W2],[f1,f2,f3,f4])

err_orig = (l,v) -> norm(K*v-l*M*v+1im*sqrt(l)*W1*v+1im*sqrt(l-σ^2)*W2*v);
err_measure = (l,v) -> err_orig(λ0+α*l,v);
λ,_,err=tiar(nep,Neig=100,displaylevel=1,maxit=100,tol=eps()*100,check_error_every=1)#,errmeasure=err_measure)

m,p=size(err);
# sort error
for j=1:p
    err[1:m,j]=sort(err[1:m,j];rev=true);
end

for j=1:p
    semilogy(1:m,err[1:m,j],color="black",linestyle="-");
end
ylim(ymax=1)
