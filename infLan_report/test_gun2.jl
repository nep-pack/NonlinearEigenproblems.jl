using NonlinearEigenproblems, Random, SparseArrays, Revise, PyPlot
import ..NEPSolver.ilan;
import ..NEPSolver.tiar;

include("../src/method_ilan.jl");
include("../src/method_tiar.jl");
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
f3 = l-> 1im*sqrt(a1*l+b1*one(l)+0im*l);
f4 = l-> 1im*sqrt(a2*l+b2*one(l)+0im*l);

# define the nep as SPMF_NEP
nep=SPMF_NEP([K-λ0*M,-α*M,W1,W2],[f1,f2,f3,f4])

# Compute the derivatives in a numerical robust way

# this function computes the first k derivatives of sqrt(ax+b) in μ
function fD_sqrt(mu,k,a,b)
    fD=zeros(ComplexF64,k)
    fD[1]=sqrt(a*mu+b)
    for j=1:k-1
        N=a*(2*j-3)
        D=2*j*(a*mu+b)
        fD[j+1]=-(N/D)*fD[j]
    end

    # manually scale with factorial since the following does not work
    #fD=broadcast(*,fD,factorial.(BigInt.(0:k-1)))
    for j=2:k for i=1:j-1
            fD[j]=fD[j]*i
    end end
    return fD
end


# mm is the number of iterations we will do
mm=50
DD=zeros(ComplexF64,2*mm,4);
DD[1,1]=1; DD[2,2]=1;
DD[:,3]=1im*fD_sqrt(0,2*mm,a1,b1);
DD[:,4]=1im*fD_sqrt(0,2*mm,a2,b2);
Dnep=DerSPMF(nep,0,DD)

err_orig = (l,v) -> norm(K*v-l*M*v+1im*sqrt(l)*W1*v+1im*sqrt(l-σ^2)*W2*v)/((norm(v))*(nK-abs(l)*nM+abs(sqrt(l))*nW1+abs(sqrt(l-σ^2))*nW2));
err_measure = (l,v) -> err_orig(λ0+α*l,v);

λ,_,err=tiar(Dnep,Neig=100,displaylevel=1,maxit=mm,tol=eps()*100,check_error_every=1,errmeasure=err_measure)

m,p=size(err);
# sort error
for j=1:p
    err[1:m,j]=sort(err[1:m,j];rev=true);
end

for j=1:p
    semilogy(1:m,err[1:m,j],color="black",linestyle="-");
end
ylim(ymax=1)
