using NonlinearEigenproblems, Random, SparseArrays, Revise, PyPlot, DelimitedFiles
import ..NEPSolver.ilan;
import ..NEPSolver.tiar;
import ..NEPSolver.iar;


include("../src/method_ilan.jl");
include("../src/method_tiar.jl");
include("../src/method_iar.jl");

nep = nep_gallery("nlevp_native_gun")

# manually loading the matrices
Av=get_Av(nep)
K=Av[1]; M=-Av[2]; W1=im*Av[3]; W2=im*Av[4]
nK=opnorm(K,1); nM=opnorm(M,1); nW1=opnorm(W1,1); nW2=opnorm(W2,1);


# define the functions
σ=108.8774; α=(300^2-200^2)/10;  λ0=250^2; # scale and shift
a1 =α; b1=λ0; a2=α; b2=λ0-σ^2
f1 = l-> one(l);
f2 = l-> l;
f3 = l-> sqrt(a1*l+b1*one(l)+0im*l);
f4 = l-> sqrt(a2*l+b2*one(l)+0im*l);

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
mm=100
DD=zeros(2*mm+2,4);
DD[1,1]=1; DD[2,2]=1;
DD[:,3]=fD_sqrt(0,2*mm+2,a1,b1);
DD[:,4]=fD_sqrt(0,2*mm+2,a2,b2);
Dnep=DerSPMF(nep,0,DD)

err_orig = (l,v) -> norm(K*v-l*M*v+sqrt(l)*W1*v+sqrt(l-σ^2)*W2*v)/((norm(v))*(nK-abs(l)*nM+abs(sqrt(l))*nW1+abs(sqrt(l-σ^2))*nW2));
err_measure = (l,v) -> err_orig(λ0+α*l,v);

λ,_,err_iar=tiar(Dnep;Neig=200,displaylevel=1,maxit=100,tol=1e-5,check_error_every=1,errmeasure=err_measure)

# plot the spectrum
λ1=λ0.+α*λ;
figure()
# region of interest
θ=range(0,stop=2*π,length=100)
plot(c.+r*cos.(θ),r*sin.(θ))
#plot(real(λ0),imag(λ0))
plot(real(λ1),imag(λ1),marker="o",markerfacecolor=:none,c=:black,linestyle=:none)       # original eigenvalues

##################################### eig approx 40 it ##############################3
mm=60
λ,_,err_iar=iar(Dnep;Neig=200,displaylevel=1,maxit=mm,tol=1e-12,check_error_every=1,errmeasure=err_measure)

# plot the spectrum
λ2=λ0.+α*λ;
plot(real(λ2),imag(λ2),marker="s",markerfacecolor=:none,c=:green,linestyle=:none)       # original eigenvalues

# run ilan
V,H,ω,HH=ilan(Dnep,Neig=80,displaylevel=1,maxit=mm,tol=1e-18,check_error_every=1,errmeasure=err_measure)
V,_,_=svd(V)



# Create a projected NEP
mm=size(V,2)
pnep=create_proj_NEP(nep,mm); # maxsize=mm
set_projectmatrices!(pnep,V,V);

n=size(K,1)
err_lifted=(λ,z)->err_measure(λ,V*z);
λ,_,err=iar(pnep;Neig=200,displaylevel=1,maxit=150,tol=1e-12,check_error_every=1,errmeasure=err_lifted)
λ3=λ0.+α*λ;
plot(real(λ3),imag(λ3),marker="*",markerfacecolor=:none,c=:red,linestyle=:none)       # original eigenvalues


# save in a file the eigenvalues and Ritz pairs
writedlm("gun_l1.csv",[real(λ1) imag(λ1)],",")
writedlm("gun_l2.csv",[real(λ2) imag(λ2)],",")
writedlm("gun_l3.csv",[real(λ3) imag(λ3)],",")
