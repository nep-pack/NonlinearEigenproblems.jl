# mm is the number of iterations we want to run
mm=100

# manually loading the matrices (since we want to shift and scale the problem)
Av=get_Av(nep);  K=Av[1]; M=-Av[2]; W1=im*Av[3]; W2=im*Av[4]
nK=opnorm(K,1); nM=opnorm(M,1); nW1=opnorm(W1,1); nW2=opnorm(W2,1);
rel_err = (l,v) -> norm(K*v-l*M*v+sqrt(l)*W1*v+sqrt(l-σ^2)*W2*v)/((norm(v))*(nK-abs(l)*nM+abs(sqrt(l))*nW1+abs(sqrt(l-σ^2))*nW2));
err_measure = (l,v) -> rel_err(λ0+α*l,v);

# define the functions (already shifted and scaled)
σ=108.8774; α=(300^2-200^2)/10;  λ0=250^2; # scale and shift paramenters
a1 =α; b1=λ0; a2=α; b2=λ0-σ^2
f1 = l-> one(l); f2 = l-> l; f3 = l-> sqrt(a1*l+b1*one(l)+0im*l); f4 = l-> sqrt(a2*l+b2*one(l)+0im*l);

# define the nep as SPMF_NEP
nep=SPMF_NEP([K-λ0*M,-α*M,W1,W2],[f1,f2,f3,f4])

# Compute the derivatives in a robust way

# this function computes the first k derivatives of sqrt(ax+b) in μ
function fD_sqrt(mu,k,a,b)
    fD=zeros(ComplexF64,k);
    fD[1]=sqrt(a*mu+b)
    for j=1:k-1
        N=a*(2*j-3); D=2*j*(a*mu+b); fD[j+1]=-(N/D)*fD[j]
    end
    for j=2:k for i=1:j-1
            fD[j]=fD[j]*i
    end end
    return fD
end

# setting up matrices containing the derivatives
DD=zeros(2*mm+2,4); DD[1,1]=1; DD[2,2]=1; DD[:,3]=fD_sqrt(0,2*mm+2,a1,b1); DD[:,4]=fD_sqrt(0,2*mm+2,a2,b2);
# create a nep defined by its derivatives in zero
Dnep=DerSPMF(nep,0,DD)
